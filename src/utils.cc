#include "utils.hpp"

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
