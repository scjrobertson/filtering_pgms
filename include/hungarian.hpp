#ifndef HUNGARIAN_HPP
#define HUNGARIAN_HPP

#include "system_constants.hpp"

double hungarianCost(Matrix<double> perf);

void hungarianStepOne(Matrix<double> & pCond,
		unsigned *stepNumber);

void hungarianStepTwo(Matrix<double> & pCond,
		Matrix<unsigned> & mask,
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		unsigned *stepNumber);

void hungarianStepThree(Matrix<unsigned> & mask, 
		ColVector<unsigned> & columnsCovered,
		unsigned *stepNumber);

void hungarianStepFour(Matrix<double> & pCond, 
		Matrix<unsigned> & mask, 
		ColVector<unsigned> & rowsCovered, 
		ColVector<unsigned> & columnsCovered, 
		unsigned *stepNumber);

Matrix<unsigned> stepFive(Matrix<unsigned> mask, 
		std::vector<unsigned> zRows, 
		std::vector<unsigned> zCols, 
		ColVector<unsigned> rowsCovered, 
		ColVector<unsigned> columnsCovered,
		unsigned *stepNumber);

Matrix<double> stepSix(Matrix<double> pCond,
		ColVector<unsigned> rowsCovered,
		ColVector<unsigned> columnsCovered,
		unsigned *stepNumber);

#endif // HUNGARIAN_HPP

