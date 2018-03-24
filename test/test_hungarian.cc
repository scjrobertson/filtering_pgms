/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for the loop assocation tests
 *************************************************************************/
#include <iostream>
#include <algorithm>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "hungarian.hpp"

class HungarianTest : public testing::Test {
	protected:
		virtual void SetUp() {
			numberOfRows = 5;
			numberOfColumns = 6;
			maxSize = std::max(numberOfColumns, numberOfRows);

			perf = gLinear::zeros<double>(numberOfRows, numberOfColumns);
			perf(0, 0) = 1.3091; perf(0, 1) = 0.2830; perf(0, 2) = 1.2613; perf(0, 3) = 0.4200; perf(0, 4) = 0.2431; perf(0, 5) = 0.7981;
			perf(1, 0) = 1.3090; perf(1, 1) = 0.4333; perf(1, 2) = 1.5169; perf(1, 3) = 0.4876; perf(1, 4) = 0.2776; perf(1, 5) = 0.4358;
			perf(2, 0) = 0.0361; perf(2, 1) = 0.5252; perf(2, 2) = 0.1509; perf(2, 3) = 0.2506; perf(2, 4) = 0.6430; perf(2, 5) = 0.9113;
			perf(3, 0) = 0.7501; perf(3, 1) = 0.7720; perf(3, 2) = 1.3429; perf(3, 3) = 0.4938; perf(3, 4) = 1.2149; perf(3, 5) = 0.5640;
			perf(4, 0) = 0.6821; perf(4, 1) = 0.2729; perf(4, 2) = 0.9993; perf(4, 3) = 0.1550; perf(4, 4) = 0.7439; perf(4, 5) = 0.5427;


			double maxValue = -1; // Assumes a distance metric is used
			for (unsigned i = 0; i < numberOfRows; i++) {
				for (unsigned j = 0; j < numberOfColumns; j++) {
					maxValue = std::max(maxValue, perf(i, j));
				} // for
			} // for

			// Create PCOnd
			pCond = Matrix<double>(maxSize, maxSize); pCond.assignToAll(maxValue);
			for (unsigned i = 0; i < numberOfRows; i++) {
				for (unsigned j = 0; j < numberOfColumns; j++) {
					pCond(i, j) = perf(i, j);
				} // for
			} // for
		}

		virtual void TearDown() {}

	protected:
		Matrix<double> perf, pCond;
		unsigned numberOfRows, numberOfColumns, maxSize;
		double maxValue;
}; // HungarianTest

TEST_F (HungarianTest, TestStepOne) {
	unsigned stepNumber = 0;
	hungarianStepOne(pCond, stepNumber);
} // TestStepOne


TEST_F (HungarianTest, TestStepTwo) {	
	unsigned stepNumber = 1;

	ColVector<unsigned> rowsCovered = ColVector<unsigned>();
	ColVector<unsigned> columnsCovered = ColVector<unsigned>();
	Matrix<unsigned> mask = Matrix<unsigned>();

	hungarianStepOne(pCond, stepNumber);
	hungarianStepTwo(pCond, mask, rowsCovered, columnsCovered, stepNumber);
} // TestStepTwo()


TEST_F (HungarianTest, TestStepThree) {	
	unsigned stepNumber = 1;

	ColVector<unsigned> rowsCovered = ColVector<unsigned>();
	ColVector<unsigned> columnsCovered = ColVector<unsigned>();
	Matrix<unsigned> mask = Matrix<unsigned>();

	hungarianStepOne(pCond, stepNumber);
	hungarianStepTwo(pCond, mask, rowsCovered, columnsCovered, stepNumber);
	hungarianStepThree(mask, columnsCovered, stepNumber);
} // TestStepThree()

TEST_F (HungarianTest, TestStepFour) {	
	unsigned stepNumber = 1;

	ColVector<unsigned> rowsCovered = ColVector<unsigned>();
	ColVector<unsigned> columnsCovered = ColVector<unsigned>();
	Matrix<unsigned> mask = Matrix<unsigned>();
	std::map<unsigned, unsigned> zRow, zCol;

	hungarianStepOne(pCond, stepNumber);
	hungarianStepTwo(pCond, mask, rowsCovered, columnsCovered, stepNumber);
	hungarianStepThree(mask, columnsCovered, stepNumber);
	hungarianStepFour(pCond, mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);
} // TestStepFour()

TEST_F (HungarianTest, TestStepSix) {	
	unsigned stepNumber = 1;

	ColVector<unsigned> rowsCovered = ColVector<unsigned>();
	ColVector<unsigned> columnsCovered = ColVector<unsigned>();
	Matrix<unsigned> mask = Matrix<unsigned>();
	std::map<unsigned, unsigned> zRow, zCol;

	hungarianStepOne(pCond, stepNumber);
	hungarianStepTwo(pCond, mask, rowsCovered, columnsCovered, stepNumber);
	hungarianStepThree(mask, columnsCovered, stepNumber);
	hungarianStepFour(pCond, mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);	
	hungarianStepSix(pCond, rowsCovered, columnsCovered, stepNumber);
} // TestStepSix()


TEST_F (HungarianTest, TestStepFourAgain) {	
	unsigned stepNumber = 1;

	ColVector<unsigned> rowsCovered = ColVector<unsigned>();
	ColVector<unsigned> columnsCovered = ColVector<unsigned>();
	Matrix<unsigned> mask = Matrix<unsigned>();
	std::map<unsigned, unsigned> zRow, zCol;

	hungarianStepOne(pCond, stepNumber);
	hungarianStepTwo(pCond, mask, rowsCovered, columnsCovered, stepNumber);
	hungarianStepThree(mask, columnsCovered, stepNumber);
	hungarianStepFour(pCond, mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);	
	hungarianStepSix(pCond, rowsCovered, columnsCovered, stepNumber);
	hungarianStepFour(pCond, mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);	
} // TestStepFourAgain()

TEST_F (HungarianTest, TestStepFive) {	
	unsigned stepNumber = 1;

	ColVector<unsigned> rowsCovered = ColVector<unsigned>();
	ColVector<unsigned> columnsCovered = ColVector<unsigned>();
	Matrix<unsigned> mask = Matrix<unsigned>();
	std::map<unsigned, unsigned> zRow, zCol;

	hungarianStepOne(pCond, stepNumber);
	hungarianStepTwo(pCond, mask, rowsCovered, columnsCovered, stepNumber);
	hungarianStepThree(mask, columnsCovered, stepNumber);
	hungarianStepFour(pCond, mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);	
	hungarianStepSix(pCond, rowsCovered, columnsCovered, stepNumber);
	hungarianStepFour(pCond, mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);	
	hungarianStepFive(mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);
} // TestStepFive()


TEST_F (HungarianTest, TestAlgorithm) {	
	double assignmentCost = hungarianCost(perf);

	std::cout << assignmentCost << std::endl;
} // TestAlgorithm()
