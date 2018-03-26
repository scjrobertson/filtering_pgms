#include <math.h>
#include <algorithm>
#include "ospa.hpp"
#include "hungarian.hpp"
#include "hellinger.hpp"


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
		for (unsigned j = 0; j < lhsSize; j++) {
			double distance = gaussianMixtureHellingerDistance(rhs[i]->w, rhs[i]->mu, rhs[i]->S, lhs[j]->w, lhs[j]->mu, lhs[j]->S);
			distanceMatrix(i, j) = pow(std::min(distance, cParameter), pParameter);
		} // for
	} // for

	double minimumCost = hungarianCost(distanceMatrix);

	unsigned maxSize = std::max(rhsSize, lhsSize);
	ospaComponents[0] = pow((1.0/maxSize)*( pow(cParameter, pParameter)*abs(rhsSize-lhsSize) + minimumCost), 1.0/pParameter );
	ospaComponents[1] = pow((1.0/maxSize)*(minimumCost), 1.0/pParameter );
	ospaComponents[0] = pow((1.0/maxSize)*(pow(cParameter, pParameter)*abs(rhsSize-lhsSize)), 1.0/pParameter );

	return ospaComponents;
} // calculateOspa()
