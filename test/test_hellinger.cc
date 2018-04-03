/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for the Hellinger distance
 *************************************************************************/
#include <iostream>
#include <algorithm>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "system_constants.hpp"
#include "hellinger.hpp"

class HellingerTest : public testing::Test {
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
		}

		virtual void TearDown() {}

	protected:
		rcptr<filters::gmm> f1, f2;
}; // HellingerTest

TEST_F (HellingerTest, TestGaussianHellingerDistance) {
	double distance = gaussianHellingerDistance(f1->mu[0], f1->S[0], f2->mu[2], f2->S[2]);
} // TestGaussianHellinger()

TEST_F (HellingerTest, TestGaussianMixtureHellingerDistance) {
	double distance = gaussianMixtureHellingerDistance(f1->w, f1->mu, f1->S, f2->w, f2->mu, f2->S);
} // TestGaussianMixtureHellinger()

TEST_F (HellingerTest, TestSingle) {

	// Mixand 1
	std::vector<double> w1 = {1};
	std::vector<ColVector<double>> mu1(1);
	std::vector<Matrix<double>> S1(1);

	mu1[0] = ColVector<double>(4);
	mu1[0][0] = 5.800314e+01; mu1[0][1] = 5.797315e+01; mu1[0][2] = 2.059062e+00; mu1[0][3] = 1.957207e+00;

	S1[0] = gLinear::zeros<double>(4, 4);
	S1[0](0, 0) = 2.275074e-01; S1[0](0, 1) = -2.104731e-04; S1[0](0, 2) = 1.477753e-01; S1[0](0, 3) = 7.957950e-05;
	S1[0](1, 0) = -2.104731e-04; S1[0](1, 1) = 2.275219e-01; S1[0](1, 2) = 7.960541e-05; S1[0](1, 3) = 1.477958e-01;
	S1[0](2, 0) = 1.477753e-01; S1[0](2, 1) = 7.960541e-05; S1[0](2, 2) = 3.103863e-01; S1[0](2, 3) = -4.859647e-05;
	S1[0](3, 0) = 7.957950e-05; S1[0](3, 1) = 1.477958e-01; S1[0](3, 2) = -4.859647e-05; S1[0](3, 3) = 3.103932e-01;

	// Mixand 2
	std::vector<double> w2 = {1};
	std::vector<ColVector<double>> mu2(1);
	std::vector<Matrix<double>> S2(1);

	mu2[0] = ColVector<double>(4);
	mu2[0][0] = 5.800355e+01; mu2[0][1] = 5.797279e+01; mu2[0][2] = 2.059906e+00; mu2[0][3] = 1.956423e+00;

	S2[0] = gLinear::zeros<double>(4, 4);
	S2[0](0, 0) = 2.270215e-01; S2[0](0, 1) = 0.000000e+00; S2[0](0, 2) = 1.479647e-01; S2[0](0, 3) = 0.000000e+00;
	S2[0](1, 0) = 0.000000e+00; S2[0](1, 1) = 2.270215e-01; S2[0](1, 2) = 0.000000e+00; S2[0](1, 3) = 1.479647e-01;
	S2[0](2, 0) = 1.479647e-01; S2[0](2, 1) = 0.000000e+00; S2[0](2, 2) = 3.102883e-01; S2[0](2, 3) = 0.000000e+00;
	S2[0](3, 0) = 0.000000e+00; S2[0](3, 1) = 1.479647e-01; S2[0](3, 2) = 0.000000e+00; S2[0](3, 3) = 3.102883e-01;

	double distance = gaussianHellingerDistance(mu1[0], S1[0], mu2[0], S2[0]);
	double distanceGM = gaussianMixtureHellingerDistance(w1, mu1, S1, w2, mu2, S2);

	std::cout << "Actual distance: " << distance << std::endl;
	std::cout << "Approximate distance: " << distanceGM << std::endl;

} // TestSingle()

TEST_F (HellingerTest, TestUHellingerForSingles) {
	std::vector<double> w1(1); w1[0] = f1->w[0];
	std::vector<ColVector<double>> mu1(1); mu1[0] = 1.0*f1->mu[0];
	std::vector<Matrix<double>> S1(1); S1[0] = 1.0*f1->S[0];

	std::vector<double> w2(1); w2[0] = f2->w[0];
	std::vector<ColVector<double>> mu2(1); mu2[0] = 1.0*f2->mu[0];
	std::vector<Matrix<double>> S2(1); S2[0] = 1.0*f2->S[0];

	/*
	double distance = gaussianHellingerDistance(f1->mu[0], f1->S[0], f2->mu[0], f2->S[0]);
	double distanceGM = gaussianMixtureHellingerDistance(w1, mu1, S1, w2, mu2, S2);

	std::cout << "Actual distance: " << distance << std::endl;
	std::cout << "Approximate distance: " << distanceGM << std::endl;
	*/
} // TestGaussianMixtureHellinger()
