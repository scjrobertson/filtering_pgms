#include <math.h>
#include <algorithm>
#include "ospa.hpp"
#include "hungarian.hpp"
#include "hellinger.hpp"

std::vector<ColVector<double>> calculateOspa(rcptr<LinearModel> model,
		std::vector<std::vector<rcptr<filters::gmm>>> & lhs,
		std::vector<std::vector<rcptr<filters::gmm>>>  & rhs) {


	unsigned rhsSize = rhs.size();
	unsigned lhsSize = lhs.size();
	unsigned ospaLength = std::max(lhs.size(), rhs.size());
	std::vector<ColVector<double>> ospa(ospaLength);

	for (unsigned i = lhsSize; i < ospaLength; i++) lhs.push_back( std::vector<rcptr<filters::gmm>>() );
	for (unsigned i = rhsSize; i < ospaLength; i++) rhs.push_back( std::vector<rcptr<filters::gmm>>() );
	for (unsigned i = 0; i < ospaLength; i++) ospa[i] = calculateOspa(lhs[i], rhs[i], model->ospaC, model->ospaP);


	return ospa;
} //  calculateOspa

ColVector<double> calculateOspa(std::vector<rcptr<filters::gmm>> lhs,
		std::vector<rcptr<filters::gmm>> rhs, 
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
			//ColVector<double> difference = rhs[i]->mu[0] - lhs[i]->mu[0];
			//double distance = pow(difference.transpose()*difference, 0.5);
			double distance = gaussianHellingerDistance(rhs[i]->mu[0], rhs[i]->S[0], lhs[j]->mu[0], lhs[j]->S[0]);
			//double distance = gaussianMixtureHellingerDistance(rhs[i]->w, rhs[i]->mu, rhs[i]->S, lhs[j]->w, lhs[j]->mu, lhs[j]->S);
			distanceMatrix(i, j) = pow(std::min(distance, cParameter), pParameter);
		} // for
	} // for

	double minimumCost = hungarianCost(distanceMatrix);

	unsigned maxSize = std::max(rhsSize, lhsSize);
	ospaComponents[0] = pow((1.0/maxSize)*( pow(cParameter, pParameter)*abs(rhsSize-lhsSize) + minimumCost), 1.0/pParameter );
	ospaComponents[1] = pow((1.0/maxSize)*(minimumCost), 1.0/pParameter );
	ospaComponents[2] = pow((1.0/maxSize)*(pow(cParameter, pParameter)*abs(rhsSize-lhsSize)), 1.0/pParameter );

	return ospaComponents;
} // calculateOspa()
