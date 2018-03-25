#include <math.h>
#include <algorithm>
#include "ospa.hpp"
#include "hungarian.hpp"


ColVector<double> calculateOspa(std::vector<rcptr<filters::gmm>> rhs,
		std::vector<rcptr<filters::gmm>> lhs, 
		double cParameter,
		double pParameter) {
	ColVector<double> ospaComponents = ColVector<double>(3); ospaComponents.assignToAll(0.0);

	unsigned rhsSize = rhs.size();
	unsigned lhsSize = lhs.size();

	// Consider all trivial cases
	if (rhsSize == 0 && lhsSize == 0) return ospaComponents;
	if (rhsSize == 0 || lhsSize == 0) {
		ospaComponents[0] = cParameter; ospaComponents[2] = cParameter;
		return ospaComponents;
	}

	// Determine indices
	Matrix<double> distanceMatrix = gLinear::zeros<double>(rhsSize, lhsSize);

	for (unsigned i = 0; i < rhsSize; i++) {
		unsigned numberOfRhsComponents = (rhs[i]->w).size();
		for (unsigned j = 0; j < lhsSize; j++) {
			unsigned numberOfLhsComponents = (lhs[i]->w).size();
			double distance = cParameter;

			if (numberOfLhsComponents == 1 && numberOfRhsComponents == 1) {
				distance = 0.0; //gaussianHellingerDistance(rhs[i]->mu[0], lhs[j]->mu[0], rhs[i]->S[0], lhs[j]->S[0]);
			} // if

			distanceMatrix(i, j) = pow(std::min(distance, cParameter), pParameter);
		} // for
	} // for

	double minimumCost = hungarianCost(distanceMatrix);


	return ospaComponents;
} // calculateOspa()
