#include <math.h>

#include "hellinger.hpp"
#include "emdw.hpp"
#include "sigmapoints.hpp"

double gaussianHellingerDistance (ColVector<double> muOne,
		Matrix<double> SOne,
		ColVector<double> muTwo,
		Matrix<double> STwo) {

	unsigned dimension = SOne.cols();
	ColVector<double> mu = muOne - muTwo;
	Matrix<double> S = SOne + STwo;

	double detP; int fail;
	Matrix<double> P = inv(S, detP, fail);

	double detSOne = det(SOne, fail);
	double detSTwo = det(STwo, fail);

	double delta = sqrt( (sqrt(detSOne*detSTwo)*detP)/pow(0.5, dimension) );
	double epsilon = exp( -0.25*(mu.transpose())*P*(mu)   );

	double hellinger = fabs(sqrt(1 - delta*epsilon));

	if (std::isnan(hellinger)) hellinger = 0.0;

	return hellinger;
} // gaussianHellingerDistance()

double gaussianMixtureHellingerDistance(std::vector<double> wOne,
		std::vector<ColVector<double>> muOne,
		std::vector<Matrix<double>> SOne,
		std::vector<double> wTwo,
		std::vector<ColVector<double>> muTwo,
		std::vector<Matrix<double>> STwo) {

	unsigned sizeOfMixtureOne = wOne.size();
	unsigned sizeOfMixtureTwo = wTwo.size();

	unsigned numberOfSigmaPoints = 2*muOne[0].size() + 1;
	double centroidWeight = 1.0/(numberOfSigmaPoints);

	std::vector<std::vector<ColVector<double>>> sigmaPointsOne(sizeOfMixtureOne);
	std::vector<std::vector<double>> sigmaPointsWeightsOne(sizeOfMixtureOne);

	std::vector<std::vector<ColVector<double>>> sigmaPointsTwo(sizeOfMixtureTwo);
	std::vector<std::vector<double>> sigmaPointsWeightsTwo(sizeOfMixtureTwo);

	// Normalise the mixtures, determine covariance matrices' inverses and determinants and get sigma points
	std::vector<double> normWOne(sizeOfMixtureOne);
	std::vector<double> normWTwo(sizeOfMixtureTwo);

	std::vector<double> detPOne(sizeOfMixtureOne);
	std::vector<double> detPTwo(sizeOfMixtureTwo);

	std::vector<Matrix<double>> POne(sizeOfMixtureOne);
	std::vector<Matrix<double>> PTwo(sizeOfMixtureTwo);

	double normConstantOne = 0.0;
	for (unsigned i = 0; i < sizeOfMixtureOne; i++) {
		// Normalising constant
		normConstantOne += wOne[i];
		// Determinant and inverse
		int fail;
		POne[i] = inv(SOne[i], detPOne[i], fail); // Could also probably do Cholesky to speed things up?
		//std::cout << "Det1: " << 1.0/detPOne[i] << std::endl;
		// Sigma Points
		cov2sp_2Dp1(muOne[i], SOne[i], centroidWeight, sigmaPointsOne[i], sigmaPointsWeightsOne[i]);
	} // for
	for (unsigned i = 0; i < sizeOfMixtureOne; i++) normWOne[i] = wOne[i]/normConstantOne;

	double normConstantTwo = 0.0;
	for (unsigned i = 0; i < sizeOfMixtureTwo; i++) {
		// Normalising constant
		normConstantTwo += wTwo[i];
		// Determinant and inverse
		int fail;
		PTwo[i] = inv(STwo[i], detPTwo[i], fail);
		//std::cout << "Det2: " << 1.0/detPTwo[i] << std::endl;
		// Sigma Points
		cov2sp_2Dp1(muTwo[i], STwo[i], centroidWeight, sigmaPointsTwo[i], sigmaPointsWeightsTwo[i]);
	} 
	for (unsigned i = 0; i < sizeOfMixtureTwo; i++) normWTwo[i] = wTwo[i]/normConstantTwo;

	// Compute the Hellinger distance
	double hSquared = 0.0;

	// Contribution from mixture one
	for (unsigned i = 0; i < sizeOfMixtureOne; i++) {
		double componentContribution = 0.0;
		for (unsigned j = 0; j < numberOfSigmaPoints; j++) {
			double mixtureOnePotential = determineMixturePotential(sigmaPointsOne[i][j], normWOne, muOne, POne, detPOne);
			double mixtureTwoPotential = determineMixturePotential(sigmaPointsOne[i][j], normWTwo, muTwo, PTwo, detPTwo);
			double rFunction = pow(sqrt(mixtureOnePotential) - sqrt(mixtureTwoPotential), 2.0)/(0.5*( mixtureOnePotential + mixtureTwoPotential));
			
			componentContribution += sigmaPointsWeightsOne[i][j]*rFunction;
			//componentContribution += (1.0/(numberOfSigmaPoints))*rFunction;
		} // for
		hSquared += 0.5*normWOne[i]*componentContribution;
	} // for

	// Contribution from mixture one
	for (unsigned i = 0; i < sizeOfMixtureTwo; i++) {
		double componentContribution = 0.0;
		for (unsigned j = 0; j < numberOfSigmaPoints; j++) {
			double mixtureOnePotential = determineMixturePotential(sigmaPointsTwo[i][j], normWOne, muOne, POne, detPOne);
			double mixtureTwoPotential = determineMixturePotential(sigmaPointsTwo[i][j], normWTwo, muTwo, PTwo, detPTwo);
			double rFunction = pow(sqrt(mixtureOnePotential) - sqrt(mixtureTwoPotential), 2.0)/(0.5*(mixtureOnePotential + mixtureTwoPotential) );

			componentContribution += sigmaPointsWeightsTwo[i][j]*rFunction;
			//componentContribution += (1.0/(numberOfSigmaPoints))*rFunction;
		} // for
		hSquared += 0.5*normWTwo[i]*componentContribution;
	} // for

	double hellingerDistance = sqrt((fabs(0.5*hSquared)));

	return hellingerDistance;
} // gaussianMixtureHellingerDistance()

double determineMixturePotential(ColVector<double> point,
		std::vector<double> w,
		std::vector<ColVector<double>> mu,
		std::vector<Matrix<double>> P,
		std::vector<double> detP) {

	double potential = 0.0;
	unsigned numberOfMixtureComponents = w.size();
	unsigned dimension = mu[0].size();

	for (unsigned i = 0; i < numberOfMixtureComponents; i++) {
		double norm = w[i]*( 1.0/( pow(detP[i], -0.5)*pow(2*M_PI, 0.5*dimension) ));
		
		ColVector<double> difference = point - mu[i];
		potential += norm*exp( -0.5*(difference.transpose())*P[i]*(difference) );
	} // for

	return potential;
} // evaluateMixtureAtPoint()
