#include "model_declaration.hpp"
#include <math.h>
#include <limits>

#include "vconstruct.hpp"
#include "mean_cov.hpp"
#include "posdef.hpp"
#include "stlvecs.hpp"
#include "vector.hpp"
#include "stdfun.hpp"

LinearModel::LinearModel() {
	// State and measurement space dimensions
	xDimension = 4;
	zDimension = 2;

	// SamplingPeriod
	samplingPeriod = 1;

	// Motion model
	A = gLinear::zeros<double>(xDimension, xDimension);
	A(0, 0) = 1; A(1, 1) = 1;
	A(0, 2) = samplingPeriod; A(1, 3) = samplingPeriod;
	A(2, 2) = 1; A(3, 3) = 1;

	u = ColVector<double>(xDimension); u.assignToAll(0.0);
	
	r0 = 0.3;
	R = gLinear::zeros<double>(xDimension, xDimension);
	R(0, 0) = pow(samplingPeriod,3)/3; R(1, 1) = pow(samplingPeriod,3)/3;
	R(0, 2) = pow(samplingPeriod,2)/2; R(1, 3) = pow(samplingPeriod,2)/2;
	R(2, 0) = pow(samplingPeriod,2)/2; R(3, 1) = pow(samplingPeriod,2)/2;
	R(2, 2) = samplingPeriod; R(3, 3) = samplingPeriod;
	R *= r0;

	// Simulation length
	simulationLength = 500;

	// Birth locations
	birthTimes = {0, 10};
	deathTimes = {100, 150};
	targetPriors.resize(2);

	targetPriors[0] = uniqptr<filters::gaussian>(new filters::gaussian);
	targetPriors[0]->id = 0;
	targetPriors[0]->w = {1.0};
	
	targetPriors[0]->mu = ColVector<double>(xDimension);
	targetPriors[0]->mu[1] = 5; targetPriors[0]->mu[2] = 5;
	targetPriors[0]->mu[2] = 1; targetPriors[0]->mu[2] = 1;

	targetPriors[0]->S = 1*R;

	targetPriors[1] = uniqptr<filters::gaussian>(new filters::gaussian);
	targetPriors[1]->id = 1;
	targetPriors[1]->w = {1.0};
	
	targetPriors[1]->mu = ColVector<double>(xDimension);
	targetPriors[1]->mu[1] = 5; targetPriors[1]->mu[2] = 5;
	targetPriors[1]->mu[2] = 1; targetPriors[1]->mu[2] = 1;

	targetPriors[1]->S = 1*R;
	
	// Meaurement model
	C = gLinear::zeros<double>(zDimension, xDimension);
	C(0, 0) = 1; C(1, 1) = 1;

	q0 = 0.3;
	Q = gLinear::zeros<double>(zDimension, zDimension);
	Q(0, 0) = 1; Q(1, 1) = 1;
	Q *= q0;

	detectionProbability = 0.95;

	// Clutter model
	observationSpaceRange.resize(zDimension);
	
	observationSpaceRange[0] = ColVector<double>(2);
	observationSpaceRange[0][0] = -50; observationSpaceRange[0][1] = 50;
	
	observationSpaceRange[1] = ColVector<double>(2);
	observationSpaceRange[1][0] = -50; observationSpaceRange[1][1] = 50;

	lambda = 60;

	// Gaussian mixture pruning parameters
	gmmComponentWeightThreshold = 1e-15;
	gmmComponentUnionDistance = std::numeric_limits<double>::infinity();
	maximumNumberOfGmmComponents = 250;

	// OSPA parameters
	ospaP = 2;
	ospaC = 5;

	this->generateGroundTruth();
} // Constructor()

void LinearModel::generateGroundTruth() {
	this->beliefs.clear();
	beliefs.resize(this->simulationLength);
	this->measurements.resize(this->simulationLength);
	this->cardinality.resize(this->simulationLength);
	unsigned N = this->birthTimes.size();

        double seed = 100;
	std::default_random_engine generator(seed);
	std::poisson_distribution<unsigned> poisson(this->lambda);
	std::uniform_real_distribution<double> uniform(0, 1.0);

	unsigned targetCounter = 0;
	// Add all targets starting at t = 0
	for (unsigned i = 0; i < N; i++) {
		if (this->birthTimes[i] == 0) {
			beliefs[0].push_back( this->targetPriors[i] );
			targetCounter++;
		} // if
	} // for

	// Declare an identity matrix
	Matrix<double> identity = gLinear::zeros<double>(xDimension, xDimension);
	for (unsigned i = 0 ; i < xDimension; i++) identity(i, i) = 1;

	for (unsigned i = 1; i < this->simulationLength; i++) {
		(this->beliefs[i]).clear();
		unsigned currentTargetNumber = beliefs[i-1].size();
		// Generate belief for exsiting targets;
		for (unsigned j = 0; j < currentTargetNumber; j++) {
			if ( i == (this->deathTimes[ beliefs[i-1][j]->id ]) ) continue;

			ColVector<double> muPred = (this->A)*beliefs[i-1][j]->mu + this->u;
			Matrix<double> covPred = (this->A)*(beliefs[i-1][j]->S)*((this->A).transpose()) + this->R;

			// Generate a measurement
			ColVector<double> zMeasurement = (this->C)*muPred + randomVector(this->zDimension, generator, 0, 4);
			if ( uniform(generator) <= this->detectionProbability ) this->measurements[i].push_back(1.0*zMeasurement);

			int fail; double detCov;
			Matrix<double> K = (covPred)*( (this->C).transpose())*inv( (this->C)*(covPred)*((this->C).transpose()) + this->Q, detCov, fail);
			ColVector<double> mu = muPred + K*( zMeasurement - (this->C)*muPred );
			Matrix<double> S = ( identity - K*(this->C) )*covPred;

			// If the target is still alive, add it to the ground truth.
			rcptr<filters::gaussian> posterior = uniqptr<filters::gaussian>(new filters::gaussian);
			posterior->id = beliefs[i-1][j]->id;
			posterior->w = 1;
			posterior->mu = 1.0*mu;
			posterior->S = 1.0*S;
			(this->beliefs[i]).push_back(posterior);
		} // for

		// Add in clutter!
		unsigned numberOfClutterMeasurements = (unsigned) poisson(generator);
		for (unsigned j = 0; j < numberOfClutterMeasurements; j++) {
			ColVector<double> ranges = ColVector<double>(2);
			ranges[0] = ((this->observationSpaceRange)[0][1] - (this->observationSpaceRange)[0][0])*uniform(generator) 
				+ (this->observationSpaceRange)[0][0];
			ranges[1] = ((this->observationSpaceRange)[1][1] - (this->observationSpaceRange)[1][0])*uniform(generator) 
				+  (this->observationSpaceRange)[1][0];

			this->measurements[i].push_back(1.0*ranges);
		} // for
		
		// Add in new targets!
		for (unsigned j = targetCounter; j < N; j++) {
			if ((i+1) == this->birthTimes[j]) {
				(this->beliefs[i]).push_back( this->targetPriors[j] );
				targetCounter++;
			} // if
		} // for

		// Add in cardinality distribubtion
		cardinality[i] = (this->beliefs[i]).size();
	} // for
} // generateGroundTruth()

std::vector<std::vector<rcptr<filters::gaussian>>> LinearModel::getGroundTruthBeliefs() const {
	return beliefs;	
} // getGroundTruthBeliefs()

std::vector<std::vector<ColVector<double>>> LinearModel::getMeasurements() const {
	return measurements;
} // getMeasurements()

std::vector<rcptr<filters::gmm>> LinearModel::getPriors(unsigned timeStep) const {
	unsigned N = targetPriors.size();
	std::vector<rcptr<filters::gmm>> newTargetPriors(N);

	for (unsigned i = 0; i < N; i++) {
		if ( birthTimes[i] == timeStep ) {
			rcptr<filters::gmm> prior = uniqptr<filters::gmm>(new filters::gmm);
			prior->w = {targetPriors[i]->w};
			prior->mu = {1.0*targetPriors[i]->mu};
			prior->S = {1.0*targetPriors[i]->S};
			newTargetPriors.push_back(prior);
		} // if
	} // for
	return newTargetPriors;
} // getPriors()

ColVector<double> randomVector(int sizeX, std::default_random_engine generator, double mean, double var){
    ColVector<double> randomVector(sizeX);
    std::normal_distribution<double> dist(mean, sqrt(var));

    for (int i=0; i < sizeX; i++) {
      randomVector[i] = dist(generator);
    }
    return randomVector;
} // normVector()


 
