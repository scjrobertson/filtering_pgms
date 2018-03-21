#ifndef OSPA_HPP
#define OSPA_HPP

#include "emdw.hpp"
#include "system_constants.hpp"
#include "model_declaration.hpp"

ColVector<double> calculateOspa(std::vector<rcptr<filters::gmm>> rhs, 
		std::vector<rcptr<filters::gmm>> lhs, 
		double c,
		double p);

double gaussianHellingerDistance (ColVector<double> muOne,
		ColVector<double> muTwo,
		Matrix<double> SOne,
		Matrix<double> STwo);

#endif // OSPA_HPP
