#ifndef HUNGARIAN_HPP
#define HUNGARIAN_HPP

#include "system_constants.hpp"

double hungarianCost(Matrix<double> perf);

Matrix<double> stepOne(Matrix<double> pCond);

Matrix<unsigned> stepTwo(Matrix<double> pCond);

ColVector<unsigned> stepThree(Matrix<unsigned> mask, 
		unsigned *stepNumber);

Matrix<unsigned> stepFour(Matrix<double> pCond, 
		Matrix<unsigned> mask, 
		ColVector<unsigned> rowsCovered, 
		ColVector<unsigned> columnsCovered, 
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

