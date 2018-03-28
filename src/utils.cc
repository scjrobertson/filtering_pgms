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
	//std::cout << "mu: " << mu << std::endl;

	// Compute the matched covariance
	for (unsigned i = 0; i < numberOfComponents; i++) {
		ColVector<double> difference = (gmm->mu[i]) - mu;
		S += normalisedWeight[i]*( (gmm->S[i]) + (difference)*(difference.transpose()) );
	} // for
	//std::cout << "S: " << S << std::endl;

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
			prunedW.push_back(gmm->w[i]);
			prunedMu.push_back(1.0*gmm->mu[i]);
			prunedS.push_back(1.0*gmm->S[i]);
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
		ColVector<double> mu = ColVector<double>(dimension);
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

void outputResults(std::vector<std::vector<ColVector<double>>> groundTruth,
		std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates,
		std::vector<ColVector<double>> ospa) {
	
	
	unsigned dimension = groundTruth[0][0].size() - 1;

	// Write out ground truth
	std::ofstream groundTruthFile;
	groundTruthFile.open("simpleGroundTruth.ini");

	unsigned numberOfTrajectories = groundTruth.size();
	for (unsigned i = 0; i < numberOfTrajectories; i++) {
		unsigned trajectoryLength = groundTruth[i].size();
		
		groundTruthFile << "[TRAJECTORY " << i+1 << "]\n";
		for (unsigned j = 0; j < trajectoryLength; j++) {
			for (unsigned k = 0; k < dimension+1; k++) groundTruthFile << groundTruth[i][j][k] << ",";
			groundTruthFile << "\n";
		} // for
		groundTruthFile << "\n";
	} // for
	groundTruthFile.close();

	// Write out state estimates
	std::ofstream stateEstimateFile;
	stateEstimateFile.open("simpleStateEstimates.ini");

	stateEstimateFile << "[DIMENSION]\n";
	stateEstimateFile << "dimension =" << dimension << "\n";

	unsigned simulationLength = stateEstimates.size();
	stateEstimateFile << "[STATE ESTIMATES]\n";
	for (unsigned i = 0; i < simulationLength; i++) {
		unsigned targetNumber = stateEstimates[i].size();
		for (unsigned j = 0; j < targetNumber; j++) {
			stateEstimateFile << i << ", ";
			unsigned numberOfGmmComponents = stateEstimates[i][j]->w.size();
			for (unsigned k = 0; k < numberOfGmmComponents; k++) {
				for (unsigned l = 0; l < dimension; l++) stateEstimateFile << stateEstimates[i][j]->mu[k][l] << ", ";
			} // for
			stateEstimateFile << "\n";
		} // for
	} // for
	stateEstimateFile.close();

	// Write out state estimates
	std::ofstream ospaFile;
	ospaFile.open("simpleOspaResults.ini");

	ospaFile << "[OSPA RESULTS]\n";
	for (unsigned i = 0; i < simulationLength; i++) {
		ospaFile << i << ", ";
		for (unsigned j = 0; j < dimension; j++) ospaFile << ospa[i][j] << ", ";
		ospaFile << "\n";
	} // for
	ospaFile.close();

} // outputResults
