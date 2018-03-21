#include "algorithmic_steps.hpp"
#include <math.h>

#include <map>
#include "factor.hpp"
#include "discretetable.hpp"
#include "clustergraph.hpp"
#include "lbp_cg.hpp"
#include "lbu_cg.hpp"

rcptr<FactorOperator> defaultInplaceNormalizer = uniqptr<FactorOperator>(new DiscreteTable_InplaceNormalize<unsigned short>);
rcptr<FactorOperator> defaultNormalizer = uniqptr<FactorOperator>(new DiscreteTable_Normalize<unsigned short>);
rcptr<FactorOperator> defaultMarginalizer = uniqptr<FactorOperator>(new DiscreteTable_Marginalize<unsigned short>);

std::vector<std::vector<rcptr<filters::gmm>>> runLinearGaussianFilter(rcptr<LinearModel> model) {
	std::vector<std::vector<rcptr<filters::gmm>>> targets; targets.clear(); targets.resize(model->simulationLength);
	std::vector<std::vector<ColVector<double>>> measurements = model->getMeasurements();

	// Add target prior for t=0
	std::vector<rcptr<filters::gmm>> targetPriors =  model->getPriors(0);	
	unsigned numberOfNewTargets = targetPriors.size();
	for (unsigned j = 0; j < numberOfNewTargets; j++) targets[0].push_back(targetPriors[j]);

	for (unsigned i = 1; i < model->simulationLength; i++) {
		//std::cout << "\nTime-step " << i << "." << std::endl;

		// Prediction
		std::vector<rcptr<filters::gmm>> predictedStates = predictMultipleTargetsLinear(model, targets[i-1]);	
		std::vector<rcptr<filters::updateComponents>> kalmanComponents = createMultipleUpdateComponentsLinear(model, predictedStates);

		// Data Association
		std::vector< std::vector<rcptr<filters::gmm>>> updateOptions = createUpdateOptionsLinear(model, predictedStates, 
				kalmanComponents, measurements[i]);
		Matrix<double> associationMatrix = createAssociationMatrix(model, measurements[i].size(), updateOptions);
		Matrix<double> updatedAssociations = loopyBeliefUpdatePropagation(model, associationMatrix);

		// Measurement Updates
		targets[i] = updateTargetStatesLinear(updateOptions, associationMatrix, updatedAssociations);

		// Add in new targets
		targetPriors =  model->getPriors(i); numberOfNewTargets = targetPriors.size();
		for (unsigned j = 0; j < numberOfNewTargets; j++) targets[i].push_back(targetPriors[j]);
	} // for

	return targets;
} // runLinearGaussianFilter()

std::vector<rcptr<filters::gmm>> predictMultipleTargetsLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> targets) {
	unsigned numberOfTargets = targets.size();
	std::vector<rcptr<filters::gmm>> predictedTargets(numberOfTargets);

	for (unsigned i = 0; i < numberOfTargets; i++) {
		unsigned numberOfComponents = (targets[i]->w).size();
		rcptr<filters::gmm> predictedGmm = uniqptr<filters::gmm>(new filters::gmm);
		predictedGmm->id = targets[i]->id;

		(predictedGmm->w).resize(numberOfComponents);
		(predictedGmm->mu).resize(numberOfComponents);
		(predictedGmm->S).resize(numberOfComponents);

		// Forward each component through the motion model
		for (unsigned j = 0; j < numberOfComponents; j++) {
			predictedGmm->w[j] = targets[i]->w[j];
			predictedGmm->mu[j] = (model->A)*(targets[i]->mu[j]) + model->u;
			predictedGmm->S[j] = (model->A)*(targets[i]->S[j])*((model->A).transpose()) + model->R;
		} // for

		predictedTargets[i] = predictedGmm;
	} // for

	return predictedTargets;
} // predictMultipleTargetsLinear()

std::vector<rcptr<filters::updateComponents>> createMultipleUpdateComponentsLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> predictedStates) {
	unsigned numberOfTargets = predictedStates.size();
	std::vector<rcptr<filters::updateComponents>> updateComponents(numberOfTargets);

	Matrix<double> eye = gLinear::zeros<double>(model->xDimension, model->xDimension);
	for (unsigned i = 0; i < model->xDimension; i++) eye(i, i) = 1.0;

	for (unsigned i = 0; i < numberOfTargets; i++) {
		unsigned numberOfComponents = (predictedStates[i]->w).size();
		rcptr<filters::updateComponents> kalmanComponents = uniqptr<filters::updateComponents>(new filters::updateComponents);

		(kalmanComponents->P).resize(numberOfComponents);
		(kalmanComponents->K).resize(numberOfComponents);
		(kalmanComponents->z).resize(numberOfComponents);
		
		(kalmanComponents->w).resize(numberOfComponents);
		(kalmanComponents->mu).resize(numberOfComponents);
		(kalmanComponents->S).resize(numberOfComponents);

		// Create an update components for each mixture component
		for (unsigned j = 0; j < numberOfComponents; j++) {
			Matrix<double> S = (model->C)*(predictedStates[i]->S[j])*((model->C).transpose()) + model->Q;

			double det; int fail;
			kalmanComponents->P[j] = inv(S, det, fail);
			kalmanComponents->K[j] = (predictedStates[i]->S[j])*((model->C).transpose())*(kalmanComponents->P[j]);
			kalmanComponents->z[j] = (model->C)*(predictedStates[i]->mu[j]);

			Matrix<double> innovation = (eye - (kalmanComponents->K[j])*(model->C));

			kalmanComponents->w[j] = (predictedStates[i]->w[j])*( 1.0/( pow(det, -0.5)*pow(2*M_PI, 0.5*model->zDimension) ));
			kalmanComponents->mu[j] = innovation*(predictedStates[i]->mu[j]);
			kalmanComponents->S[j] = innovation*(predictedStates[i]->S[j]);
		} // for

		updateComponents[i] = kalmanComponents;
	} // for

	return updateComponents;
} // predictMultipleTargetsLinear()

std::vector< std::vector<rcptr<filters::gmm>>> createUpdateOptionsLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> predictedStates,
		std::vector<rcptr<filters::updateComponents>> kalmanComponents,
		std::vector<ColVector<double>> z
		) {
	
	unsigned numberOfTargets = kalmanComponents.size();
	unsigned numberOfMeasurements = z.size();
	std::vector<std::vector<rcptr<filters::gmm>>> updateOptions(numberOfTargets);

	for (unsigned i = 0; i < numberOfTargets; i++) {
		unsigned numberOfMixtureComponents = (predictedStates[i]->w).size();
		updateOptions[i].resize(numberOfMeasurements+1);

		// Add in missed measurement option
		rcptr<filters::gmm> predictedComponent = uniqptr<filters::gmm>(new filters::gmm);
		(predictedComponent->w).resize(numberOfMixtureComponents); // Lots of this, put a method in the struct?
		(predictedComponent->mu).resize(numberOfMixtureComponents);
		(predictedComponent->S).resize(numberOfMixtureComponents);

		for (unsigned j = 0; j < numberOfMixtureComponents; j++) {
			predictedComponent->w[j] = (1 - model->detectionProbability)*(predictedStates[i]->w[j])/(model->observationSpaceVolume);
			predictedComponent->mu[j] = predictedStates[i]->mu[j];
			predictedComponent->S[j] = predictedStates[i]->S[j];
		} // for
		updateOptions[i][0] = predictedComponent;

		// Add in the rest of the options
		for (unsigned j = 0; j < numberOfMeasurements; j++) {
			rcptr<filters::gmm> updatedComponent = uniqptr<filters::gmm>(new filters::gmm);
			(updatedComponent->w).resize(numberOfMixtureComponents);
			(updatedComponent->mu).resize(numberOfMixtureComponents);
			(updatedComponent->S).resize(numberOfMixtureComponents);
			
			for (unsigned k = 0; k < numberOfMixtureComponents; k++) {
				ColVector<double> difference = z[j] - kalmanComponents[i]->z[k];
				double likelihood = exp(-0.5*(difference.transpose()*(kalmanComponents[i]->P[k])*difference));

				updatedComponent->w[k] = (model->detectionProbability)*(kalmanComponents[i]->w[k])*(likelihood);
				updatedComponent->mu[k] = kalmanComponents[i]->mu[k] + kalmanComponents[i]->K[k]*z[j];
				updatedComponent->S[k] = kalmanComponents[i]->S[k];
			} // for
			updateOptions[i][j+1] = updatedComponent;
		} // for
	} // for
	return updateOptions;
} // updateMultipleLinear()

Matrix<double> createAssociationMatrix(rcptr<LinearModel> linearModel, 
		unsigned numberOfMeasurements,
		std::vector< std::vector<rcptr<filters::gmm>>> updateComponents) {
	unsigned numberOfTargets = updateComponents.size();
	Matrix<double> associationMatrix = gLinear::zeros<double>(numberOfTargets, numberOfMeasurements + 1);
	
	for (unsigned i = 0; i < numberOfTargets; i++) {
		unsigned numberOfMixtureComponents = (updateComponents[i][0]->w).size();
		for (unsigned j = 0; j < numberOfMeasurements+1; j++) {
			double likelihood = 0;
			for (unsigned k = 0; k < numberOfMixtureComponents; k++) likelihood += updateComponents[i][j]->w[k];
			associationMatrix[i][j] = likelihood;
		} // for
	} // for

	return associationMatrix;
} // createAssociationMatrix()

Matrix<double> loopyBeliefUpdatePropagation (rcptr<LinearModel> linearModel,
		Matrix<double> associationMatrix) {
	unsigned numberOfTargets = associationMatrix.rows();
	unsigned numberOfUpdateOptions = associationMatrix.cols();
	
	Matrix<double> gatedAssociations = gLinear::zeros<double>(numberOfTargets, numberOfUpdateOptions);
	std::vector<std::vector<unsigned short>> domains(numberOfTargets);
	
	Matrix<double> updatedAssociations = gLinear::zeros<double>(numberOfTargets, numberOfUpdateOptions);
	
	std::vector<rcptr<Factor>> marginalAssociations(numberOfTargets);
	std::vector<rcptr<Factor>> associationPriors(numberOfTargets);

	// Gate the associations
	for (unsigned i = 0; i < numberOfTargets; i++) {
		double normalisingConstant = 0;
		for (unsigned j = 0; j < numberOfUpdateOptions; j++) {
			if (associationMatrix[i][j] >= 0) {
				gatedAssociations[i][j] = associationMatrix[i][j];
				normalisingConstant += associationMatrix[i][j];
				domains[i].push_back( (unsigned short) j);
			} // if
		} // for

		// Normalise the association probabilities
		for (unsigned j = 0; j < numberOfUpdateOptions; j++) gatedAssociations[i][j] /= normalisingConstant;

		// Create the marginal association factors
		std::map<std::vector<unsigned short>, FProb> sparseProbsMarginal;
		std::map<std::vector<unsigned short>, FProb> sparseProbsPrior;
		
		rcptr<std::vector<unsigned short>> marginalDom = uniqptr<std::vector<unsigned short>>( new std::vector<unsigned short> (domains[i]));
		rcptr<std::vector<unsigned short>> priorDom = uniqptr<std::vector<unsigned short>>( new std::vector<unsigned short> (domains[i]));

		for (unsigned j = 0; j < domains[i].size(); j++) {
			sparseProbsMarginal[{domains[i][j]}] = gatedAssociations[i][j]; 
			sparseProbsPrior[{domains[i][j]}] = 1; 
		} // for

		marginalAssociations[i] = uniqptr<Factor> (new DiscreteTable<unsigned short>(
					emdw::RVIds{i}, {marginalDom},
					0.0, sparseProbsMarginal, 
					0.0, 0.0,
					false, defaultMarginalizer,
					defaultInplaceNormalizer, defaultNormalizer) );

		associationPriors[i] = uniqptr<Factor> (new DiscreteTable<unsigned short>(
					emdw::RVIds{i}, {priorDom},
					0.0, sparseProbsPrior, 
					0.0, 0.0,
					false, defaultMarginalizer,
					defaultInplaceNormalizer, defaultNormalizer) );
	} // for


	// Check for disjoint clusters
	std::vector<bool> isNotDisjoint(numberOfTargets);
	for (unsigned i = 0; i < numberOfTargets; i++) isNotDisjoint[i] = false;

	std::vector<rcptr<Factor>> clusters;
	for (unsigned i = 0; i < numberOfTargets; i++) {
		unsigned short aSize = domains[i].size();

		for (unsigned j = i+1; j < numberOfTargets; j++) {
			if (haveIntersectingDomains(domains[i], domains[j])) {
				isNotDisjoint[i] = true; isNotDisjoint[j] = true;

				rcptr<Factor> product = associationPriors[i]->absorb(associationPriors[j]);
				rcptr<DiscreteTable<unsigned short>> cancel = std::dynamic_pointer_cast<DiscreteTable<unsigned short>>(product);
				emdw::RVIds scope = {i, j};

				for (unsigned short k = 1; k < aSize; k++) cancel->setEntry(scope, emdw::RVVals{ k, k }, 0 );
				product->inplaceNormalize();

				clusters.push_back( product );
			} // if
		} // for
	} // for

	// Add in disjoint clusters
	for (unsigned i = 0; i < numberOfTargets; i++) {
		if ( isNotDisjoint[i] == false ) clusters.push_back( uniqptr<Factor> (marginalAssociations[i]->copy()) );
	} // for

	// Absorb in marginal associations
	unsigned numberOfClusters = clusters.size();

	for (unsigned i = 0; i < numberOfTargets; i++) {
		
		emdw::RVIds assocationVariable = marginalAssociations[i]->getVars();

		for (unsigned j = 0; j < numberOfClusters; j++) {
			emdw::RVIds vars = clusters[j]->getVars();

			if ( vars.size() == 1 ) {
				if ( vars[0] == i ) {
					clusters[j]->inplaceAbsorb(marginalAssociations[i]);
					clusters[j]->inplaceNormalize();
				} // if
			} else if (vars.size() == 2) {
				if ( (vars[0] == i || vars[1] == i) ) {
					clusters[j]->inplaceAbsorb(marginalAssociations[i]);
					clusters[j]->inplaceNormalize();
				} // if
			} // if
		} // for
	} // for


	// Create the cluster graph
	rcptr<ClusterGraph> clusterGraph = uniqptr<ClusterGraph>(new ClusterGraph(clusters));
	std::map<Idx2, rcptr<Factor>> msgs; msgs.clear();
	MessageQueue msgQ; msgQ.clear();

	//
	std::vector<rcptr<Factor>> finalAssociationProbs(numberOfTargets);
	try {
		// Pass messages until convergence
		unsigned nMsg = loopyBU_CG(*clusterGraph, msgs, msgQ, 0.33);

		for (unsigned i = 0; i < numberOfTargets; i++) {
			emdw::RVIds nodeVars = {i};
			finalAssociationProbs[i] = queryLBU_CG(*clusterGraph, msgs, nodeVars )->normalize();
		} // for
	} // try() 
	catch (const char* msg) {
		std::cerr << msg << std::endl;
		throw;
	} // catch()
	catch (const std::string& msg) {
		std::cerr << msg << std::endl;
		throw;
	} // catch()
	catch (const std::exception& e) {
		std::cerr << "Unhandled exception: " << e.what() << std::endl;
		throw e;
	} // catch()
	catch(...) {
		std::cerr << "An unknown exception / error occurred\n";
		throw;
	} // catch()

	// Create updated association matrix
	for (unsigned i = 0; i < numberOfTargets; i++) {
		unsigned domainSize = domains[i].size();
		double normalisingConstant = 0;
		rcptr<DiscreteTable<unsigned short>> marginal = std::dynamic_pointer_cast<DiscreteTable<unsigned short>>(finalAssociationProbs[i]);

		for (unsigned j = 0; j < domainSize; j++) {
			double potential = marginal->potentialAt(emdw::RVIds{i}, emdw::RVVals{ domains[i][j] });
			if (potential > 0.0) {
				updatedAssociations[i][j] = potential;
				normalisingConstant += potential;
			} // if
		} // for
		
		for (unsigned j = 0; j < numberOfUpdateOptions; j++) updatedAssociations[i][j] /= normalisingConstant;
	} // for

	return updatedAssociations;
} // loopyBeliefUpdatePropagation()

std::vector<rcptr<filters::gmm>> updateTargetStatesLinear(std::vector<std::vector<rcptr<filters::gmm>>> updateOptions,
		Matrix<double> associationMatrix,
		Matrix<double> updatedAssociations) {

	unsigned numberOfTargets = associationMatrix.rows();
	unsigned numberOfUpdateOptions = associationMatrix.cols();

	std::vector<rcptr<filters::gmm>> updatedTargets(numberOfTargets);

	for (unsigned i = 0; i < numberOfTargets; i++) {
		rcptr<filters::gmm> gmm = uniqptr<filters::gmm>(new filters::gmm);
		(gmm->w).clear(); (gmm->mu).clear(); (gmm->S).clear();

		unsigned numberOfMixtureComponents = (updateOptions[i][0]->w).size();
		for (unsigned j = 0; j < numberOfUpdateOptions; j++) {
			if (updatedAssociations[i][j] > 0.0) { // Assumes loopyBeliefUpdatePropagation gates final probabilities
				double associationWeight = updatedAssociations[i][j]/associationMatrix[i][j];

				for (unsigned k = 0; k < numberOfMixtureComponents; k++) {
					(gmm->w).push_back(associationWeight*(updateOptions[i][j]->w[k]));
					(gmm->mu).push_back(1.0*updateOptions[i][j]->mu[k]);
					(gmm->S).push_back(1.0*updateOptions[i][j]->S[k]);
				} // for
			} // if
		} // for
		updatedTargets[i] = weakMarginalisation(gmm);
	} // for

	return updatedTargets;
} // updatedTargetStatesLinear()
