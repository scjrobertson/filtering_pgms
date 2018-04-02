#ifndef SYSTEM_CONSTANTS_HPP
#define SYSTEM_CONSTANTS_HPP

#include <vector>
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "factor.hpp"
#include "discretetable.hpp"

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

	struct cfm {
		unsigned id = 0;
		std::vector<double> g;
		std::vector<ColVector<double>> h;
		std::vector<Matrix<double>> K;
	};

} // namespace

#endif // SYSTEM_CONSTANTS_HPP
