#ifndef MODEL_DECLARATION_HPP
#define MODEL_DECLARATION_HPP

#include <vector>
#include <random>
#include <ctime>
#include "system_constants.hpp"
#include "genvec.hpp"
#include "genmat.hpp"

ColVector<double> randomVector(int sizeX, std::default_random_engine generator, double mean, double variance);

class LinearModel {

	public:
		LinearModel();

	public:
		unsigned xDimension;
		unsigned zDimension;
	
	public:
		double samplingPeriod;

	public:
		Matrix<double> A;
		ColVector<double> u;

		double r0;
		Matrix<double> R;

	public:
		unsigned simulationLength;
		std::vector<unsigned> birthTimes;
		std::vector<unsigned> deathTimes;
		std::vector<rcptr<filters::gmm>> targetPriors;

	public:
		Matrix<double> C;

		double q0;
		Matrix<double> Q;

		double detectionProbability;

	public:
		std::vector<ColVector<double>> observationSpaceRange;
		int lambda;
		double observationSpaceVolume;

	public:
		std::vector<std::vector<ColVector<double>>> groundTruth;
		std::vector<std::vector<rcptr<filters::gmm>>> beliefs;
		std::vector<unsigned> cardinality;

	public:
		std::vector<std::vector<ColVector<double>>> measurements;

	private:
		void generateGroundTruth();

	public:
		std::vector<std::vector<ColVector<double>>> getIndividualGroundTruthTrajectories() const;
		std::vector<std::vector<ColVector<double>>> getGroundTruth() const;
		std::vector<std::vector<rcptr<filters::gmm>>> getGroundTruthBeliefs() const;
		std::vector<std::vector<ColVector<double>>> getMeasurements() const;

	public:
		std::vector<rcptr<filters::gmm>> getPriors(unsigned timeStep) const;

	public:
		double gmmComponentWeightThreshold;
		double gmmComponentUnionDistance;
		int maximumNumberOfGmmComponents;

	public:
		double ospaP;
		double ospaC;
};

#endif // MODEL_DECLARATION_HPP
