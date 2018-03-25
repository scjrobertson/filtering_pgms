#ifndef OSPA_HPP
#define OSPA_HPP

#include "emdw.hpp"
#include "system_constants.hpp"
#include "hellinger.hpp"

ColVector<double> calculateOspa(std::vector<rcptr<filters::gmm>> rhs, 
		std::vector<rcptr<filters::gmm>> lhs, 
		double c,
		double p);
#endif // OSPA_HPP
