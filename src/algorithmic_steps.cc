#include "algorithmic_steps.hpp"
#include <math.h>

#include "factor.hpp"
#include "discretetable.hpp"
#include "clustergraph.hpp"
#include "lbp_cg.hpp"
#include "lbu_cg.hpp"

std::vector<std::vector<rcptr<filters::gmm>>> runLinearGaussianFilter(rcptr<LinearModel> model) {
	std::vector<std::vector<rcptr<filters::gmm>>> targets; targets.clear(); targets.resize(model->simulationLength);
	std::vector<std::vector<rcptr<filters::gaussian>>> groundTruth = model->getGroundTruthBeliefs();
	std::vector<std::vector<ColVector<double>>> measurements = model->getMeasurements();

	// Add target prior for t=0
	std::vector<rcptr<filters::gmm>> targetPriors =  model->getPriors(0);	
	unsigned numberOfNewTargets = targetPriors.size();
	for (unsigned j = 0; j < numberOfNewTargets; j++) targets[0].push_back(targetPriors[j]);

	for (unsigned i = 1; i < model->simulationLength; i++) {
		// Prediction
		std::vector<rcptr<filters::gmm>> predictedStates = predictMultipleTargetsLinear(model, targets[i-1]);	
		std::vector<rcptr<filters::updateComponents>> kalmanComponents = createMultipleUpdateComponentsLinear(model, predictedStates);

		// Measurement update
		std::vector< std::vector<rcptr<filters::gmm>>> updateOptions = updateMultipleLinear(model, predictedStates, kalmanComponents, measurements[i]);
		Matrix<double> associationMatrix = createAssociationMatrix(model, measurements[i].size(), updateOptions);
		Matrix<double> updatedAssociations = loopyBeliefUpdatePropagation(model, associationMatrix);

		// Re-allocate the existing targets
		unsigned currentTargetNumber = targets[i-1].size(); targets[i].resize(currentTargetNumber);
		for (unsigned j = 0; j < currentTargetNumber; j++) targets[i][j] = predictedStates[j];

		// Add in new targets
		std::vector<rcptr<filters::gmm>> targetPriors =  model->getPriors(i);	
		unsigned numberOfNewTargets = targetPriors.size();
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

			kalmanComponents->w[j] = (predictedStates[i]->w[j])*( 1.0/( pow(1.0/det, 0.5)*pow(M_PI, 0.5*model->zDimension) ));
			kalmanComponents->mu[j] = innovation*(predictedStates[i]->mu[j]);
			kalmanComponents->S[j] = innovation*(predictedStates[i]->S[j]);
		} // for

		updateComponents[i] = kalmanComponents;
	} // for

	return updateComponents;
} // predictMultipleTargetsLinear()

std::vector< std::vector<rcptr<filters::gmm>>> updateMultipleLinear(rcptr<LinearModel> model, 
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
			predictedComponent->w[j] = (1 - model->detectionProbability)*predictedStates[i]->w[j];
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
	std::vector<std::vector<unsigned>> domains(numberOfTargets);
	
	Matrix<double> updatedAssociations = gLinear::zeros<double>(numberOfTargets, numberOfUpdateOptions);

	// Gate the associations
	for (unsigned i = 0; i < numberOfTargets; i++) {
		double normalisingConstant = 0;
		for (unsigned j = 0; j < numberOfUpdateOptions; j++) {
			if (associationMatrix[i][j] >= 0) {
				gatedAssociations[i][j] = associationMatrix[i][j];
				normalisingConstant += associationMatrix[i][j];
				domains[i].push_back(j);
			} // if
		} // for
		for (unsigned j = 0; j < numberOfUpdateOptions; j++) gatedAssociations[i][j] /= normalisingConstant;
	} // for
	
	return updatedAssociations;
} // loopyBeliefUpdatePropagation()
