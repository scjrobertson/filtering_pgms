#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>

rcptr<filters::gmm> weakMarginalisation(rcptr<filters::gmm> gmm) {
	unsigned numberOfComponents = (gmm->w).size();
	unsigned dimension = (gmm->mu[0]).size();
	
	double w = 0;
	std::vector<double> normalisedWeight(numberOfComponents);
	ColVector<double> mu = ColVector<double>(dimension); mu.assignToAll(0.0);
	Matrix<double> S = gLinear::zeros<double>(dimension, dimension);

	// Determine total weight of the mixture
	for (unsigned i = 0; i < numberOfComponents; i++) w += (gmm->w[i]);

	// Compute the mached mean
	for (unsigned i = 0; i < numberOfComponents; i++) {
		normalisedWeight[i] = (gmm->w[i])/w;
		mu += normalisedWeight[i]*(gmm->mu[i]);
	} // for

	// Compute the matched covariance
	for (unsigned i = 0; i < numberOfComponents; i++) {
		ColVector<double> difference = (gmm->mu[i]) - mu;
		S += normalisedWeight[i]*( (gmm->S[i]) + (difference)*(difference.transpose()) );
	} // for

	// Allocate the weak marginal
	rcptr<filters::gmm> weakMarginal = uniqptr<filters::gmm>(new filters::gmm);
	weakMarginal->id = gmm->id;
	weakMarginal->w = {1}; // Maybe preserve the orginal mass?
	weakMarginal->mu = {mu};
	weakMarginal->S = {S};

	return weakMarginal;
} // weakMarginalisation()

rcptr<filters::gmm> gaussianMixturePruning(rcptr<filters::gmm> gmm,
		double componentWeightThreshold,
		double componentUnionDistance,
		unsigned maximumNumberOfMixtureComponents) {
	unsigned numberOfComponents = gmm->w.size();

	if (numberOfComponents == 0) return gmm;

	unsigned dimension = (gmm->mu[0]).size();

	// Remove insignificant components
	std::vector<double> prunedW; prunedW.clear();
	std::vector<ColVector<double>> prunedMu; prunedMu.clear();
	std::vector<Matrix<double>> prunedS; prunedS.clear();

	for (unsigned i = 0; i < numberOfComponents; i++) {
		if (gmm->w[i] >= componentWeightThreshold) {
			prunedW.push_back(gmm->w[i]); //std::cout << "w: " << gmm->w[i] << std::endl;
			prunedMu.push_back(1.0*gmm->mu[i]); //std::cout << "mu: " <<  gmm->mu[i] << std::endl;
			prunedS.push_back(1.0*gmm->S[i]); //std::cout << "S: " << gmm->S[i] << std::endl;
		} else {
		} // if
	} // for

	// Merge closely spaced components
	std::vector<double> mergedW; mergedW.clear();
	std::vector<ColVector<double>> mergedMu; mergedMu.clear();
	std::vector<Matrix<double>> mergedS; mergedS.clear();

	unsigned numberOfRemainingComponents = prunedW.size();
	while(numberOfRemainingComponents != 0) {
		std::vector<unsigned> mergedIndices; mergedIndices.clear();
		// Get biggest component
 		double maxWeight = -1;
		unsigned maxIndex = -1;
		for (unsigned i = 0; i < numberOfRemainingComponents; i++) maxWeight = std::max(maxWeight, prunedW[i]);
		for (unsigned i = 0; i < numberOfRemainingComponents; i++) {
			if (prunedW[i] == maxWeight) {
				maxIndex = i;
				break;
			} // if
		} // for
		// Merge closely spaced components
		double w = 0;
		ColVector<double> mu = ColVector<double>(dimension); mu.assignToAll(0.0);
		Matrix<double> S = gLinear::zeros<double>(dimension, dimension);

		double det; int fail;
		Matrix<double> P = inv(prunedS[maxIndex], det, fail);

		for (unsigned i = 0; i < numberOfRemainingComponents; i++) {
			ColVector<double> difference = prunedMu[i] - prunedMu[maxIndex];
			double mahalanobis = (difference.transpose())*P*(difference);

			if (mahalanobis <= componentUnionDistance ) {
				w += prunedW[i];
				mu += prunedW[i]*prunedMu[i];
				S += prunedW[i]*( prunedS[i] + (difference)*(difference.transpose()) );
				mergedIndices.push_back(i);
			} // if
		} // for

		// Adjust for the weights
		mergedW.push_back(w);
		mergedMu.push_back(mu/w);
		mergedS.push_back(S/w);

		// Re-allocate
		unsigned numberOfMergedComponents = mergedIndices.size();
		std::vector<double> tempW; tempW.clear();
		std::vector<ColVector<double>> tempMu; tempMu.clear();
		std::vector<Matrix<double>> tempS; tempS.clear();

		//std::cout << "mergedIndices: " << mergedIndices << "\n" << std::endl;

		for (unsigned i = 0; i < numberOfRemainingComponents; i++) {
			bool merged = false;
			for (unsigned j = 0; j < numberOfMergedComponents; j++) {
				if (i == mergedIndices[j]) {
					merged = true;
					break;
				} // if
			} // for
			if (!merged) {
				tempW.push_back(prunedW[i]);
				tempMu.push_back(1.0*prunedMu[i]);
				tempS.push_back(1.0*prunedS[i]);
			} // if
		} // for
		numberOfRemainingComponents = tempW.size();
		prunedW.clear(); prunedMu.clear(); prunedS.clear();
		prunedW.resize(numberOfRemainingComponents); prunedMu.resize(numberOfRemainingComponents); prunedS.resize(numberOfRemainingComponents);

		for (unsigned i = 0; i < numberOfRemainingComponents; i++) {
			prunedW[i] = tempW[i];
			prunedMu[i] = 1.0*tempMu[i];
			prunedS[i] = 1.0*tempS[i];
		} // for
	} // while
	
	rcptr<filters::gmm> prunedAndMergedGaussianMixture = uniqptr<filters::gmm>(new filters::gmm);

	unsigned numberOfMergedComponents = mergedW.size();
	if ( numberOfMergedComponents > maximumNumberOfMixtureComponents ) {
		// Cap the total number of mixands
		std::vector<double> tempW; tempW.clear();
		std::vector<ColVector<double>> tempMu; tempMu.clear();
		std::vector<Matrix<double>> tempS; tempS.clear();

		std::vector<size_t> sortedIndices = sortIndices(mergedW, std::greater<double>());

		for (unsigned i = 0; i < maximumNumberOfMixtureComponents; i++) {
			tempW.push_back(mergedW[sortedIndices[i]]);
			tempMu.push_back(mergedMu[sortedIndices[i]]);
			tempS.push_back(mergedS[sortedIndices[i]]);
		} // for

		// Reassign
		(prunedAndMergedGaussianMixture->w) = tempW; (prunedAndMergedGaussianMixture->mu) = tempMu; 
		(prunedAndMergedGaussianMixture->S) = tempS;
		
	} else {
		(prunedAndMergedGaussianMixture->w) = mergedW; (prunedAndMergedGaussianMixture->mu) = mergedMu; 
		(prunedAndMergedGaussianMixture->S) = mergedS;
	} // if

	return prunedAndMergedGaussianMixture;
} // gaussianMixturePruning()

bool haveIntersectingDomains(std::vector<unsigned short> a, std::vector<unsigned short> b) {
	unsigned aSize = a.size();
	unsigned bSize = b.size();
	bool doIntersect = false;

	for (unsigned i = 0; i < aSize; i++) {
		for (unsigned j = 0; j < bSize; j++) {
			if (a[i] == b[i]) {  
				doIntersect = true;
				break;
			} // if
		}
	} // for

	return doIntersect;
} // haveIntersectingDomains

Matrix<double> measurementToTargetTransform(Matrix<double> associationMatrix) {
	unsigned numberOfTargets = associationMatrix.rows();
	unsigned numberOfMeasurements = associationMatrix.cols()-1; 

	Matrix<double> transformedMatrix = gLinear::zeros<double>(numberOfMeasurements, numberOfTargets+1);
	
	// Account for measurement being clutter!
	for(unsigned i = 0; i < numberOfMeasurements; i++) transformedMatrix(i, 0) = associationMatrix(0, 0);

	// Transpose the remaining enteries
	for (unsigned i = 0; i < numberOfMeasurements; i++) {
		for (unsigned j = 1; j < numberOfTargets+1; j++) transformedMatrix(i, j) = associationMatrix(j-1, i+1);
	} // for

	return transformedMatrix;
} // measurementToTargetTransform()

rcptr<filters::cfm> convertGmmToCfm( rcptr<filters::gmm> gmm ) {
	rcptr<filters::cfm> convert = uniqptr<filters::cfm>(new filters::cfm);
	
	unsigned numberOfComponents = gmm->w.size();
	convert->g = std::vector<double>(numberOfComponents);
	convert->h = std::vector<ColVector<double>>(numberOfComponents);
	convert->K = std::vector<Matrix<double>>(numberOfComponents);

	for (unsigned i = 0; i < numberOfComponents; i++) {
		convert->g[i] = log(gmm->w[i]);
		double det; int fail;
		convert->K[i] = inv(gmm->S[i], det, fail);
		convert->h[i] = (convert->K[i])*(gmm->mu[i]);
	} // for

	return convert;
} // convertGmmToCfm()

rcptr<filters::gmm> convertCfmToGmm(rcptr<filters::cfm> cfm) {
	rcptr<filters::gmm> convert = uniqptr<filters::gmm>(new filters::gmm);
	
	unsigned numberOfComponents = cfm->g.size();
	convert->w = std::vector<double>(numberOfComponents);
	convert->mu = std::vector<ColVector<double>>(numberOfComponents);
	convert->S = std::vector<Matrix<double>>(numberOfComponents);

	for (unsigned i = 0; i < numberOfComponents; i++) {
		convert->w[i] = exp(cfm->g[i]);
		double det; int fail;
		convert->S[i] = inv(cfm->K[i], det, fail);
		convert->mu[i] = (convert->S[i])*(cfm->h[i]);
	} // for

	return convert;
} // convertGmmToCfm()

void outputResults(rcptr<LinearModel> model,
		std::vector<std::vector<ColVector<double>>> groundTruth,
		std::vector<std::vector<ColVector<double>>> measurements,
		std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates,
		std::vector<ColVector<double>> ospa,
		std::vector<unsigned> cardinality) {
	
	unsigned xDimension = model->xDimension;
	unsigned zDimension = model->zDimension;
	unsigned lengthOfSimulation = measurements.size();
	Matrix<double> C = 1.0*model->C;

	// Write out ground truth
	std::ofstream groundTruthFile;
	groundTruthFile.open("plotting/data/simpleGroundTruth.ini");

	unsigned numberOfTrajectories = groundTruth.size();

	groundTruthFile << "[SIMULATION INFO]\n";
	groundTruthFile << "xDimension= " << xDimension << "\n";
	groundTruthFile << "zDimension= " << zDimension << "\n";
	groundTruthFile << "simulationLength= " << lengthOfSimulation << "\n";
	groundTruthFile << "numberOfGroundTruthTrajectories= " << numberOfTrajectories << "\n";
	groundTruthFile << "observationModel=";
	for (unsigned i = 0; i < zDimension; i++) {
		for (unsigned j = 0; j < xDimension; j++) groundTruthFile << C(i, j) << ", ";
	} // for
	groundTruthFile << "\n";

	for (unsigned i = 0; i < numberOfTrajectories; i++) {
		unsigned trajectoryLength = groundTruth[i].size();
		
		groundTruthFile << "[TRAJECTORY " << i+1 << "]\n";
		groundTruthFile << "trajectory =";
		for (unsigned j = 0; j < trajectoryLength; j++) {
			groundTruthFile << groundTruth[i][j][0] << ", ";
			for (unsigned k = 0; k < xDimension; k++) groundTruthFile << groundTruth[i][j][k+1] << ", ";
		} // for
		groundTruthFile << "\n";
	} // for
	groundTruthFile.close();

	// Write measurements
	std::ofstream measurementsFile;
	measurementsFile.open("plotting/data/simpleMeasurements.ini");

	for (unsigned i = 0; i < lengthOfSimulation; i++) {
		measurementsFile << "[MEASUREMENTS " << i+1 << "]\n";
		measurementsFile << "measurments = ";
		unsigned numberOfMeasurements = measurements[i].size();
		
		for (unsigned j = 0; j < numberOfMeasurements; j++) {
			for (unsigned k = 0; k < zDimension; k++) measurementsFile << measurements[i][j][k] << ", ";
		} // for
		measurementsFile << "\n";
	} // for
	measurementsFile.close();

	// Write out state estimates
	std::ofstream stateEstimateFile;
	stateEstimateFile.open("plotting/data/simpleStateEstimates.ini");
	
	unsigned simulationLength = stateEstimates.size();
	for (unsigned i = 0; i < simulationLength; i++) {
		unsigned targetNumber = stateEstimates[i].size();

		stateEstimateFile << "[STATE ESTIMATES " << i+1 << "]\n";
		stateEstimateFile << "numberOfTargets = " << targetNumber << "\n";

		for (unsigned j = 0; j < targetNumber; j++) {
			stateEstimateFile << "target" << j+1 << "= ";
			unsigned numberOfGmmComponents = stateEstimates[i][j]->w.size();
			for (unsigned k = 0; k < numberOfGmmComponents; k++) {
				for (unsigned l = 0; l < xDimension; l++) stateEstimateFile << stateEstimates[i][j]->mu[k][l] << ", ";
			} // for
			stateEstimateFile << "\n";
		} // for
	} // for
	stateEstimateFile.close();

	// Write out state estimates
	std::ofstream ospaFile;
	ospaFile.open("plotting/data/simpleOspaResults.ini");

	ospaFile << "[OSPA]\n";
	ospaFile << "ospa = ";
	for (unsigned i = 0; i < simulationLength; i++) {
		for (unsigned j = 0; j < 3; j++) ospaFile << ospa[i][j] << ", ";
	} // for
	ospaFile << "\ncardinality = ";
	for (unsigned i = 0; i < simulationLength; i++) {
		ospaFile << cardinality[i] << ", ";
	} // for
	ospaFile.close();

} // outputResults

void printOspaTrials(std::vector<std::vector<std::vector<ColVector<double>>>> ospa) { 	
	unsigned numberOfTrials = ospa.size();
	if (numberOfTrials == 0) return;
	unsigned numberOfSimulationsPerTrial = ospa[0].size();
	
	std::ofstream trialFile;
	trialFile.open("plotting/data/jpdafClutterTrials.csv");

	trialFile << numberOfTrials << "\n";
	trialFile << numberOfSimulationsPerTrial << "\n";
	trialFile << 50 << "\n";
	trialFile << 1 << "\n";

	for (unsigned i = 0; i < numberOfTrials; i++) {
		for (unsigned j = 0; j < numberOfSimulationsPerTrial; j++) {
			unsigned simulationLength = ospa[i][j].size();
			for (unsigned k = 0; k < simulationLength; k++) trialFile << ospa[i][j][k][1] << ",";
			trialFile << "\n";
		} // for
	} // for
	
	trialFile.close();
} // printOspaTrials()
