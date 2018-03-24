#ifndef HUNGARIAN_HPP
#define HUNGARIAN_HPP

#include "system_constants.hpp"
#include <map>

double hungarianCost(Matrix<double> perf);

void hungarianStepOne(Matrix<double> & pCond,
		unsigned & stepNumber);

void hungarianStepTwo(Matrix<double> & pCond,
		Matrix<unsigned> & mask,
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		unsigned & stepNumber);

void hungarianStepThree(Matrix<unsigned> & mask, 
		ColVector<unsigned> & columnsCovered,
		unsigned & stepNumber);

void hungarianStepFour(Matrix<double> & pCond, 
		Matrix<unsigned> & mask, 
		ColVector<unsigned> & rowsCovered, 
		ColVector<unsigned> & columnsCovered, 
	 	std::map<unsigned, unsigned> & zRow,
		std::map<unsigned, unsigned> & zCol,
		unsigned & stepNumber);

void hungarianStepFive(Matrix<unsigned> & mask, 
		ColVector<unsigned> & rowsCovered, 
		ColVector<unsigned> & columnsCovered, 
		std::map<unsigned, unsigned> & zRow, 
		std::map<unsigned, unsigned> & zCol,
		unsigned & stepNumber);

void hungarianStepSix(Matrix<double> & pCond,
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		unsigned & stepNumber);

#endif // HUNGARIAN_HPP

