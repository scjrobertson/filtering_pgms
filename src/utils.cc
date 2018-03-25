#include "utils.hpp"
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

		//std::cout << "w[" << i << "]: " << gmm->w[i] << std::endl;
		//std::cout << "mu[" << i << "]: " << gmm->mu[i] << std::endl;
		mu += normalisedWeight[i]*(gmm->mu[i]);
	} // for
	//std::cout << "mu: " << mu << std::endl;

	// Compute the matched covariance
	for (unsigned i = 0; i < numberOfComponents; i++) {
		ColVector<double> difference = (gmm->mu[i]) - mu;
		//std::cout << "S[" << i << "]: " << gmm->S[i] << std::endl;
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
 		double maxWeight = std::numeric_limits<double>::infinity();
		unsigned maxIndex = -1;
		for (unsigned i = 0; i < numberOfRemainingComponents; i++) std::max(maxWeight, prunedW[i]);
		for (unsigned i = 0; i < numberOfRemainingComponents; i++) if (prunedW[i] == maxWeight) maxIndex = i;
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

		for (unsigned i = 0; i < numberOfRemainingComponents; i++) {
			bool merged = false;
			for (unsigned j = 0; j < numberOfMergedComponents; j++) if (i == mergedIndices[j]) merged = true;
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
} // for
