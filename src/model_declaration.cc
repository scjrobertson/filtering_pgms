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
	simulationLength = 50;

	// Birth locations
	birthTimes = {0};
	deathTimes = {50};
	targetPriors.resize(1);

	// Target 1
	targetPriors[0] = uniqptr<filters::gmm>(new filters::gmm);
	targetPriors[0]->id = 0;
	targetPriors[0]->w = {1.0};
	
	(targetPriors[0]->mu).resize(1);
	targetPriors[0]->mu[0] = ColVector<double>(xDimension);
	targetPriors[0]->mu[0][0] = -40; targetPriors[0]->mu[0][1] = 40;
	targetPriors[0]->mu[0][2] = 1.0; targetPriors[0]->mu[0][3] = -2.0;

	(targetPriors[0]->S).resize(1);
	targetPriors[0]->S[0] = gLinear::zeros<double>(xDimension, xDimension);
	targetPriors[0]->S[0](0, 0) = 1.0; targetPriors[0]->S[0](1, 1) = 1.0;
	targetPriors[0]->S[0](2, 2) = 1.0; targetPriors[0]->S[0](3, 3) = 1.0;
	
	// Target 2
	/*
	targetPriors[1] = uniqptr<filters::gmm>(new filters::gmm);
	targetPriors[1]->id = 1;
	targetPriors[1]->w = {1.0};
	
	(targetPriors[1]->mu).resize(1);
	targetPriors[1]->mu[0] = ColVector<double>(xDimension);
	targetPriors[1]->mu[0][0] = -40.0; targetPriors[1]->mu[0][1] = -40.0;
	targetPriors[1]->mu[0][2] = 2.0; targetPriors[1]->mu[0][3] = 2.0;

	(targetPriors[1]->S).resize(1);
	targetPriors[1]->S[0] = gLinear::zeros<double>(xDimension, xDimension);
	targetPriors[1]->S[0](0, 0) = 1.0; targetPriors[1]->S[0](1, 1) = 1.0;
	targetPriors[1]->S[0](2, 2) = 1.0; targetPriors[1]->S[0](3, 3) = 1.0;
	*/
	
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
	observationSpaceRange[0][0] = -100; observationSpaceRange[0][1] = 100;
	
	observationSpaceRange[1] = ColVector<double>(2);
	observationSpaceRange[1][0] = -100; observationSpaceRange[1][1] = 100;

	observationSpaceVolume = 1e4;

	lambda = 60;

	// Gaussian mixture pruning parameters
	gmmComponentWeightThreshold = 1e-30;
	gmmComponentUnionDistance = 4;//std::numeric_limits<double>::infinity();
	maximumNumberOfGmmComponents = 100;

	// OSPA parameters
	ospaP = 2;
	ospaC = 1;

	this->generateGroundTruth();
} // Constructor()

void LinearModel::generateGroundTruth() {
	this->groundTruth.clear();
	groundTruth.resize(this->simulationLength);
	this->beliefs.clear();
	beliefs.resize(this->simulationLength);
	this->measurements.resize(this->simulationLength);
	this->cardinality.resize(this->simulationLength);
	unsigned numberOfTargets = this->birthTimes.size();

	std::default_random_engine generator(time(0));
	std::poisson_distribution<unsigned> poisson(this->lambda);
	std::uniform_real_distribution<double> uniform(0, 1.0);

	unsigned targetCounter = 0;
	// Add all targets starting at t = 0
	for (unsigned i = 0; i < numberOfTargets; i++) {
		if (this->birthTimes[i] == 0) {
			(this->beliefs[0]).push_back( this->targetPriors[i] );
			if ((this->targetPriors[i]->mu).size() > 0) (this->groundTruth[0]).push_back(this->targetPriors[i]->mu[0]);
			targetCounter++;
		} // if
	} // for

	cardinality[0] = (this->beliefs[0]).size();

	// Declare an identity matrix
	Matrix<double> identity = gLinear::zeros<double>(xDimension, xDimension);
	for (unsigned i = 0 ; i < xDimension; i++) identity(i, i) = 1;

	for (unsigned i = 1; i < this->simulationLength; i++) {
		(this->beliefs[i]).clear(); (this->groundTruth[i]).clear();
		unsigned currentTargetNumber = beliefs[i-1].size();
		// Generate belief for exsiting targets;
		for (unsigned j = 0; j < currentTargetNumber; j++) {
			if ( i == (this->deathTimes[ beliefs[i-1][j]->id ]) ) continue;

			ColVector<double> muPred = (this->A)*(beliefs[i-1][j]->mu[0]) + this->u;
			Matrix<double> covPred = (this->A)*(beliefs[i-1][j]->S[0])*((this->A).transpose()) + this->R;

			// Generate a measurement
			ColVector<double> muTruth = (this->A)*groundTruth[i-1][j] + this->u;
			ColVector<double> zMeasurement = (this->C)*muTruth + randomVector(this->zDimension, generator, 0, 0.25);
			if ( uniform(generator) <= this->detectionProbability ) this->measurements[i].push_back(1.0*zMeasurement);

			int fail; double detCov;
			Matrix<double> K = (covPred)*( (this->C).transpose())*inv( (this->C)*(covPred)*((this->C).transpose()) + this->Q, detCov, fail);
			ColVector<double> mu = muPred + K*( zMeasurement - (this->C)*muPred );
			Matrix<double> S = ( identity - K*(this->C) )*covPred;

			// If the target is still alive, add it to the ground truth beliefs.
			rcptr<filters::gmm> posterior = uniqptr<filters::gmm>(new filters::gmm);
			posterior->id = beliefs[i-1][j]->id;
			posterior->w = {1};
			posterior->mu = {1.0*mu};
			posterior->S = {1.0*S};
			(this->beliefs[i]).push_back(posterior);

			// Add ground truths
			(this->groundTruth[i]).push_back(muTruth);
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
		for (unsigned j = targetCounter; j < numberOfTargets; j++) {
			if ((i+1) == this->birthTimes[j]) {
				if ((this->targetPriors[j]->mu).size() > 0) (this->groundTruth[i]).push_back(this->targetPriors[j]->mu[0]);
				(this->beliefs[i]).push_back(this->targetPriors[j]);
				targetCounter++;
			} // if
		} // for

		// Add in cardinality distribubtion
		cardinality[i] = (this->beliefs[i]).size();
	} // for
} // generateGroundTruth()


std::vector<std::vector<ColVector<double>>> LinearModel::getIndividualGroundTruthTrajectories() const {
	unsigned numberOfTrajectories = targetPriors.size();
	std::vector<std::vector<ColVector<double>>> individualTrajectories(numberOfTrajectories);

	for (unsigned i = 0; i < numberOfTrajectories; i++) {
		unsigned trajectoryLength = deathTimes[i] - birthTimes[i];
		std::vector<ColVector<double>> trajectory(trajectoryLength);
		std::vector<ColVector<double>> exportTrajectory(trajectoryLength);
		
		// Initial position and time
		exportTrajectory[0] = ColVector<double>(xDimension+1); exportTrajectory[0].assignToAll(0.0);
		exportTrajectory[0][0] = birthTimes[i];
		trajectory[0] = 1.0*targetPriors[i]->mu[0];
		for (unsigned j = 0; j < xDimension; j++) exportTrajectory[0][j+1] = trajectory[0][j];

		for (unsigned j = 1; j < trajectoryLength; j++) {
			exportTrajectory[j] = ColVector<double>(xDimension+1); exportTrajectory[j].assignToAll(0.0);
			exportTrajectory[j][0] = birthTimes[i] + j;

			trajectory[j] = A*trajectory[j-1] + u;
			for (unsigned k = 0; k < xDimension; k++) exportTrajectory[j][k+1] = trajectory[j][k];	
		} // for
		individualTrajectories[i] = exportTrajectory;
	} // for

	return individualTrajectories;	
} // getGroundTruth()


std::vector<std::vector<ColVector<double>>> LinearModel::getGroundTruth() const {
	return groundTruth;	
} // getGroundTruth()

std::vector<std::vector<rcptr<filters::gmm>>> LinearModel::getGroundTruthBeliefs() const {
	return beliefs;	
} // getGroundTruthBeliefs()

std::vector<std::vector<ColVector<double>>> LinearModel::getMeasurements() const {
	return measurements;
} // getMeasurements()

std::vector<rcptr<filters::gmm>> LinearModel::getPriors(unsigned timeStep) const {

	unsigned N = targetPriors.size();
	std::vector<rcptr<filters::gmm>> newTargetPriors; newTargetPriors.clear();

	for (unsigned i = 0; i < N; i++) {
		if ( birthTimes[i] == timeStep ) {
			if (targetPriors.size() > 0 ) {  
				rcptr<filters::gmm> prior = uniqptr<filters::gmm>(new filters::gmm);
				prior->w = {targetPriors[i]->w[0]};
				prior->mu = {1.0*targetPriors[i]->mu[0]};
				prior->S = {1.0*targetPriors[i]->S[0]};
				newTargetPriors.push_back(prior);
			}
		} // if
	} // for

	return newTargetPriors;
} // getPriors()

std::vector<unsigned> LinearModel::getCardinality() const {
	return cardinality;
} // getCardinality()

ColVector<double> randomVector(int sizeX, std::default_random_engine generator, double mean, double var){
    ColVector<double> randomVector(sizeX);
    std::normal_distribution<double> dist(mean, sqrt(var));

    for (int i=0; i < sizeX; i++) {
      randomVector[i] = dist(generator);
    }
    return randomVector;
} // normVector()


 
