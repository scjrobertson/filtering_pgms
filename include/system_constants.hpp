#ifndef SYSTEM_CONSTANTS_HPP
#define SYSTEM_CONSTANTS_HPP

#include <vector>
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"

namespace filters {

	struct gaussian {
		unsigned id = 0;
		double w = 1;
		ColVector<double> mu;
		Matrix<double> S;
	};

	struct gmm {
		unsigned id = 0;
		std::vector<double> w;
		std::vector<ColVector<double>> mu;
		std::vector<Matrix<double>> S;
	};

	struct updateComponents {
		std::vector<Matrix<double>> P;
		std::vector<Matrix<double>> K;
		std::vector<ColVector<double>> z;

		std::vector<double> w;
		std::vector<ColVector<double>> mu;
		std::vector<Matrix<double>> S;
	};
} // namespace

#endif // SYSTEM_CONSTANTS_HPP
