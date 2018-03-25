#ifndef HELLINGER_HPP
#define HELLINGER_HPP

#include "system_constants.hpp"

double gaussianHellingerDistance (ColVector<double> muOne,
		Matrix<double> SOne,
		ColVector<double> muTwo,
		Matrix<double> STwo);

double gaussianMixtureHellingerDistance (std::vector<double> wOne,
		std::vector<ColVector<double>> muOne,
		std::vector<Matrix<double>> SOne,
		std::vector<double> wTwo,
		std::vector<ColVector<double>> muTwo,
		std::vector<Matrix<double>> STwo);

double determineMixturePotential (ColVector<double> point,
		std::vector<double> w,
		std::vector<ColVector<double>> mu,
		std::vector<Matrix<double>> P,
		std::vector<double> detP);

#endif // HELLINGER_HPP
