#include "algorithmic_steps.hpp"
#include <math.h>

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
		std::vector<rcptr<filters::gmm>> predictedTargets = predictMultipleTargetsLinear(model, targets[i-1]);		

		// Re-allocate the existing targets
		unsigned currentTargetNumber = targets[i-1].size(); targets[i].resize(currentTargetNumber);
		for (unsigned j = 0; j < currentTargetNumber; j++) targets[i][j] = predictedTargets[j];

		// Add in new targets
		std::vector<rcptr<filters::gmm>> targetPriors =  model->getPriors(i);	
		unsigned numberOfNewTargets = targetPriors.size();
		for (unsigned j = 0; j < numberOfNewTargets; j++) targets[i].push_back(targetPriors[j]);

		std::cout << "t: " << i << ", new target number: " << numberOfNewTargets << std::endl;
		std::cout << "t: " << i << ", current target number: " << targets[i].size() << std::endl;
	} // for

	return targets;
} // runLinearGaussianFilter()

std::vector<rcptr<filters::gmm>> predictMultipleTargetsLinear(rcptr<LinearModel> model, std::vector<rcptr<filters::gmm>> targets) {
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

std::vector<rcptr<filters::updateComponents>> createMultipleUpdateComponentsLinear(rcptr<LinearModel> model, std::vector<rcptr<filters::gmm>> predictedStates) {
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

			kalmanComponents->w[j] = (predictedStates[i]->w[j])*(1.0/( pow(det, 0.5)*pow(M_PI, 0.5*model->zDimension) ));
			kalmanComponents->mu[j] = innovation*(predictedStates[i]->mu[j]);
			kalmanComponents->S[j] = innovation*(predictedStates[i]->S[j]);
		} // for

		updateComponents[i] = kalmanComponents;
	} // for

	return updateComponents;
} // predictMultipleTargetsLinear()
