/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for the utility functions
 *************************************************************************/
#include <iostream>
#include <algorithm>
#include <limits>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "sortindices.hpp"
#include "system_constants.hpp"
#include "utils.hpp"

class UtilsTest : public testing::Test {
	protected:
		virtual void SetUp() {
			// Mixture 1
			f1 = uniqptr<filters::gmm>(new filters::gmm);

			(f1->w).resize(5); 
			(f1->w[0]) = 0.0747; (f1->w[1]) = 0.2507;
			(f1->w[2]) = 0.1346; (f1->w[3]) = 0.1934;
			(f1->w[4]) = 0.3467;
			
			(f1->mu).resize(5); 
			for (unsigned i = 0; i < 5; i++) (f1->mu[i]) = ColVector<double>(4);
			
			(f1->mu[0][0]) = 0.3963; (f1->mu[0][1]) = 0.7051;
			(f1->mu[0][2]) = 0.5586; (f1->mu[0][3]) = 0.7566;

			(f1->mu[1][0]) = 0.9955; (f1->mu[1][1]) = 0.9624;
			(f1->mu[1][2]) = 0.5351; (f1->mu[1][3]) = 0.9639;
			
			(f1->mu[2][0]) = 0.1156; (f1->mu[2][1]) = 0.0514;
			(f1->mu[2][2]) = 0.3043; (f1->mu[2][3]) = 0.5802;

			(f1->mu[3][0]) = 0.5310; (f1->mu[3][1]) = 0.9012;
			(f1->mu[3][2]) = 0.5406; (f1->mu[3][3]) = 0.4320;

			(f1->mu[4][0]) = 0.5427; (f1->mu[4][1]) = 0.7124;
			(f1->mu[4][2]) = 0.0167; (f1->mu[4][3]) = 0.8009;
			
			(f1->S).resize(5);
			for (unsigned i = 0; i < 5; i++) (f1->S[i]) = gLinear::zeros<double>(4, 4);

			(f1->S[0])(0, 0) = 0.1696; (f1->S[0])(1, 1) = 0.2788;
			(f1->S[0])(2, 2) = 0.1982; (f1->S[0])(3, 3) = 0.1951;

			(f1->S[1])(0, 0) = 0.3268; (f1->S[1])(1, 1) = 0.8803;
			(f1->S[1])(2, 2) = 0.4711; (f1->S[1])(3, 3) = 0.4040;

			(f1->S[2])(0, 0) = 0.1792; (f1->S[2])(1, 1) = 0.9689;
			(f1->S[2])(2, 2) = 0.4075; (f1->S[2])(3, 3) = 0.8445;

			(f1->S[3])(0, 0) = 0.6153; (f1->S[3])(1, 1) = 0.3766;
			(f1->S[3])(2, 2) = 0.8772; (f1->S[3])(3, 3) = 0.7849;

			(f1->S[4])(0, 0) = 0.4650; (f1->S[4])(1, 1) = 0.8140;
			(f1->S[4])(2, 2) = 0.8984; (f1->S[4])(3, 3) = 0.4292;

			// Mixture 2
			f2 = uniqptr<filters::gmm>(new filters::gmm);

			(f2->w).resize(3); 
			(f2->w[0]) = 0.1643; (f2->w[1]) = 0.8192;
			(f2->w[2]) = 0.0165;

			(f2->mu).resize(3); 
			for (unsigned i = 0; i < 3; i++) (f2->mu[i]) = ColVector<double>(4);
			
			(f2->mu[0][0]) = 2.3343; (f2->mu[0][1]) = 2.5966;
			(f2->mu[0][2]) = 2.9020; (f2->mu[0][3]) = 2.7021;

			(f2->mu[1][0]) = 2.3775; (f2->mu[1][1]) = 2.7350;
			(f2->mu[1][2]) = 2.9541; (f2->mu[1][3]) = 2.5428;
			
			(f2->mu[2][0]) = 2.5401; (f2->mu[2][1]) = 2.3111;
			(f2->mu[2][2]) = 2.0712; (f2->mu[2][3]) = 2.1820;

			(f2->S).resize(3);
			for (unsigned i = 0; i < 3; i++) (f2->S[i]) = gLinear::zeros<double>(4, 4);

			(f2->S[0])(0, 0) = 0.9150; (f2->S[0])(1, 1) = 0.6427;
			(f2->S[0])(2, 2) = 0.0014; (f2->S[0])(3, 3) = 0.0304;

			(f2->S[1])(0, 0) = 0.2085; (f2->S[1])(1, 1) = 0.4550;
			(f2->S[1])(2, 2) = 0.1273; (f2->S[1])(3, 3) = 0.0086;

			(f2->S[2])(0, 0) = 0.7271; (f2->S[2])(1, 1) = 0.3541;
			(f2->S[2])(2, 2) = 0.7804; (f2->S[2])(3, 3) = 0.4367;

			weightThreshold = 1e-15;
			unionDistance = std::numeric_limits<double>::infinity();
			maximumNumberOfComponents = 1;
		}

		virtual void TearDown() {}

	protected:
		rcptr<filters::gmm> f1, f2;
		double weightThreshold, unionDistance;
		unsigned maximumNumberOfComponents;
}; // UtilsTest

TEST_F (UtilsTest, TestGMMPruningInfinite) {
	double componentUnionDistance = 0.5;

	rcptr<filters::gmm> prunedGaussianMixture = gaussianMixturePruning(f1, weightThreshold, componentUnionDistance, maximumNumberOfComponents);
	rcptr<filters::gmm> weakMarginal = weakMarginalisation(f1);

	std::cout << "prunedGaussianMixture: " << std::endl;
	std::cout << (prunedGaussianMixture->w[0]) << std::endl;
	std::cout << (prunedGaussianMixture->mu[0]) << std::endl;
	std::cout << (prunedGaussianMixture->S[0]) << std::endl;

	std::cout << "WeakMarginal: " << std::endl;
	std::cout << (weakMarginal->w[0]) << std::endl;
	std::cout << (weakMarginal->mu[0]) << std::endl;
	std::cout << (weakMarginal->S[0]) << std::endl;
} // TestGMMPruning()
