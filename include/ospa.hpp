#ifndef OSPA_HPP
#define OSPA_HPP

#include "emdw.hpp"
#include "system_constants.hpp"
#include "hellinger.hpp"
#include "model_declaration.hpp"

std::vector<ColVector<double>> calculateOspa(rcptr<LinearModel> model,
		std::vector<std::vector<rcptr<filters::gmm>>> & lhs,
		std::vector<std::vector<rcptr<filters::gmm>>> & rhs);

ColVector<double> calculateOspa(std::vector<rcptr<filters::gmm>> lhs, 
		std::vector<rcptr<filters::gmm>> rhs, 
		double c,
		double p);
#endif // OSPA_HPP
