/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for the loopy belief propagation
 *************************************************************************/
#include <iostream>
#include <algorithm>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "system_constants.hpp"
#include "model_declaration.hpp"
#include "algorithmic_steps.hpp"

class LBPTest : public testing::Test {
	protected:
		virtual void SetUp() {
			model = uniqptr<LinearModel>(new LinearModel(20, 0.95));
			A = gLinear::zeros<double>(4, 7);
    			A(0, 0) = 0.8147; A(0, 1) = 0.6324; A(0, 2) = 0.9575; A(0, 3) = 0.9572; A(0, 4) = 0.4218; A(0, 5) = 0.6557; A(0, 6) = 0.6787;
    			A(1, 0) = 0.9058; A(1, 1) = 0.0975; A(1, 2) = 0.9649; A(1, 3) = 0.4854; A(1, 4) = 0.9157; A(1, 5) = 0.0357; A(1, 6) = 0.7577;
   		 	A(2, 0) = 0.1270; A(2, 1) = 0.2785; A(2, 2) = 0.1576; A(2, 3) = 0.8003; A(2, 4) = 0.7922; A(2, 5) = 0.8491; A(2, 6) = 0.7431;
    			A(3, 0) = 0.9134; A(3, 1) = 0.5469; A(3, 2) = 0.9706; A(3, 3) = 0.1419; A(3, 4) = 0.9595; A(3, 5) = 0.9340; A(3, 6) = 0.3922;
		}

		virtual void TearDown() {}

	protected:
		rcptr<LinearModel> model;
		Matrix<double> A, B;
}; // HellingerTest

TEST_F (LBPTest, TestLoopyBeliefPropagation) {
	Matrix<double> D = loopyBeliefPropagation(model, A);

	std::cout << "D: " << D << std::endl;
} // TestGaussianHellinger()

