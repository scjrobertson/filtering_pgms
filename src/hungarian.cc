#include "hungarian.hpp"
#include <algorithm>
#include <limits>

double hungarianCost(Matrix<double> perf) {

	unsigned numberOfRows = perf.rows();
	unsigned numberOfColumns = perf.cols();
	unsigned maxSize = std::max(numberOfRows, numberOfColumns);

	Matrix<double> matching = gLinear::zeros<double>(numberOfRows, numberOfColumns);
	Matrix<double> edge = gLinear::zeros<double>(maxSize, maxSize);

	// Get max value
	double maxValue = -1; // Assumes a distance metric is used
	for (unsigned i = 0; i < numberOfRows; i++) {
		for (unsigned j = 0; j < numberOfColumns; j++) {
			maxValue = std::max(maxValue, perf(i, j));
		} // for
	} // for

	// Create PCOnd
	Matrix<double> pCond(maxSize, maxSize); pCond.assignToAll(maxValue);
	for (unsigned i = 0; i < numberOfRows; i++) {
		for (unsigned j = 0; j < numberOfColumns; j++) {
			pCond(i, j) = perf(i, j);
		} // for
	} // for
	

	bool exitFlag = true;
	unsigned stepNumber = 1;

	ColVector<unsigned> rowsCovered = ColVector<unsigned>();
	ColVector<unsigned> columnsCovered = ColVector<unsigned>();
	Matrix<unsigned> mask = Matrix<unsigned>();
	std::map<unsigned, unsigned> zRow, zCol;

	while(exitFlag) {
		//std::cout << "StepNumber: " << stepNumber << std::endl;
		switch(stepNumber) {
			case 1:
				hungarianStepOne(pCond, stepNumber);
				break;
			case 2:
				hungarianStepTwo(pCond, mask, rowsCovered, columnsCovered, stepNumber);
				//std::cout << mask << std::endl;
				break;
			case 3:
				hungarianStepThree(mask, columnsCovered, stepNumber);
				//std::cout << mask << std::endl;
				break;
			case 4:
				hungarianStepFour(pCond, mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);
				//std::cout << mask << std::endl;
				break;
			case 5:
				hungarianStepFive(mask, rowsCovered, columnsCovered, zRow, zCol, stepNumber);
				//std::cout << mask << std::endl;
				break;
			case 6:
				hungarianStepSix(pCond, rowsCovered, columnsCovered, stepNumber);
				break;
			case 7:
				exitFlag = false;
				break;
		} // switch
	} // while

	double assignmentCost = 0.0;
	for (unsigned i = 0; i < numberOfRows; i++) {
		for (unsigned j = 0; j < numberOfColumns; j++) {
			if (mask(i, j) == 1) assignmentCost += perf(i, j);
		} // for
	} // for

	return assignmentCost;
} // hungarianCost()

void hungarianStepOne(Matrix<double> & pCond,
		unsigned & stepNumber) {
	unsigned numberOfRows = pCond.rows();
	unsigned numberOfColumns = pCond.cols();

	for (unsigned i = 0; i < numberOfRows; i++) {
		double minValue = std::numeric_limits<double>::infinity();
		for (unsigned j = 0; j < numberOfColumns; j++) minValue = std::min(minValue, pCond(i, j));
		for (unsigned j = 0; j < numberOfColumns; j++) pCond(i, j) = pCond(i,j) - minValue;
	} // for

	stepNumber = 2;
} // stepOne()

void hungarianStepTwo(Matrix<double> & pCond,
		Matrix<unsigned> & mask,
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		unsigned & stepNumber) {
	unsigned numberOfRows = pCond.rows();
	unsigned numberOfColumns = pCond.cols();

	mask.resize(numberOfRows, numberOfColumns); mask.assignToAll(0);
	rowsCovered.resize(numberOfRows); rowsCovered.assignToAll(0);
	columnsCovered.resize(numberOfColumns); columnsCovered.assignToAll(0);

	for (unsigned i = 0; i < numberOfRows; i++) {
		for (unsigned j = 0; j < numberOfColumns; j++) {
			if ( pCond(i, j) == 0 && rowsCovered[i] == 0 && columnsCovered[j] == 0 ) {
				mask(i, j) = 1;
				rowsCovered[i] = 1;
				columnsCovered[j] = 1;
			} // if
		} // for
	} // for

	rowsCovered.assignToAll(0);
	columnsCovered.assignToAll(0);

	stepNumber = 3;
} // stepTwo()

void hungarianStepThree(Matrix<unsigned> & mask,
		ColVector<unsigned> & columnsCovered,
		unsigned & stepNumber) {
	unsigned numberOfColumns = mask.cols();
	columnsCovered.resize(numberOfColumns); columnsCovered.assignToAll(0);
	unsigned totalSum = 0.0;

	for (unsigned i = 0; i < numberOfColumns; i++) {
		for (unsigned j = 0; j < numberOfColumns; j++) {
			if (mask(i, j)) {
				columnsCovered[j]++;
				totalSum++;
			} // if
		} // for
	} // for

	if (totalSum == numberOfColumns) stepNumber = 7;
	else stepNumber = 4;
} // stepThree()

void hungarianStepFour(Matrix<double> & pCond, 
		Matrix<unsigned> & mask, 
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		std::map<unsigned, unsigned> & zRow,
		std::map<unsigned, unsigned> & zCol,
		unsigned & stepNumber) {

	unsigned pSize = pCond.cols();
	bool zFlag = true;

	zRow.clear(); zCol.clear();

	while(zFlag) {
		int row = -1; int col = -1;
		bool exitFlag = true;
		unsigned i = 0; unsigned j = 0;

		while (exitFlag) {
			if (pCond(i, j) == 0 && rowsCovered[i] == 0 && columnsCovered[j] == 0 ) {
				row = i; 
				col = j;
				exitFlag = false;
			} // if
			j++;
			if (j > pSize) {  
				j = 1; i++;
			}
			if (i > pSize) exitFlag = false;
		} // while

		if (row == -1) {
			stepNumber = 6;
			zFlag = false;
			zRow[0] = 0; 
			zCol[0] = 0;
		} else {
			mask(row, col) = 2;

			std::vector<unsigned> indices;
			unsigned indexSum = 0;
			for (unsigned i = 0; i < pSize; i++) {
				if (mask(row, i) == 1) {
					indices.push_back(i);
					indexSum += i;
				} // if
			} // for

			if (indexSum != 0) {
				rowsCovered[row] = 1;
				for (unsigned i = 0; i < indices.size(); i++) columnsCovered[indices[i]] = 0;
			} else {
				stepNumber = 5;
				zFlag = false;
				zRow[0] = row;
				zCol[0] = col;
			} // if
		} // if
	} // while
} // stepFour()

void hungarianStepFive(Matrix<unsigned> & mask,
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		std::map<unsigned, unsigned> & zRow,
		std::map<unsigned, unsigned> & zCol,
		unsigned & stepNumber) {
	
	unsigned pSize = mask.cols();

	bool zFlag = true;
	unsigned i = 0;
	while (zFlag) {
		int rowIndex = -1;
		for (unsigned j = 0; j < pSize; j++) if ( mask(j, zCol[i]) == 1 ) rowIndex = j;

		if (rowIndex > -1) {
			i++;
			zRow[i] = (unsigned) rowIndex;
			zCol[i] = zCol[i-1];
		} else {
			zFlag = false;
		} // if
		
		if (zFlag == true) {
			int columnIndex = -1;
			for (unsigned j = 0; j < pSize; j++) if ( mask(zRow[i], j) == 2 ) columnIndex = j;
			i++;
			zRow[i] = zRow[i-1];
			zCol[i] = (unsigned) columnIndex;
		} // if
	} // while

	for (unsigned i = 0; i < zRow.size(); i++) {
		if (mask(zRow[i], zCol[i]) == 1 ) mask(zRow[i], zCol[i]) = 0;
		else mask(zRow[i], zCol[i]) = 1;
	} // for

	zRow.clear(); zCol.clear();

	for (unsigned i = 0; i < pSize; i++) {
		for (unsigned j = 0; j < pSize; j++) if (mask(i, j) == 2) mask(i, j) = 0;
	}

	stepNumber = 3;
} // stepFive()

void hungarianStepSix(Matrix<double> & pCond,
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		unsigned & stepNumber) {
	
	unsigned rowsSize = rowsCovered.size();
	unsigned columnsSize = columnsCovered.size();

	std::vector<unsigned> a, b, c;
	for (unsigned i = 0; i < rowsSize; i++) {
		if(rowsCovered[i] == 0) a.push_back(i);
		if(rowsCovered[i] == 1) c.push_back(i);
	} 
	for (unsigned i = 0; i < columnsSize; i++) if(columnsCovered[i] == 0) b.push_back(i);

	unsigned pSize = pCond.cols();
	unsigned aSize = a.size(); 
	unsigned bSize = b.size();
	unsigned cSize = c.size();

	double minVal = std::numeric_limits<double>::infinity();
	for (unsigned i = 0; i < aSize; i++) {
		for (unsigned j = 0; j < bSize; j++) {
			minVal = std::min(minVal, pCond(a[i], b[j]));
		} // for
	} // for

	for (unsigned i = 0; i < cSize; i++ ) {
		for (unsigned j = 0; j < pSize; j++) pCond(c[i], j) += minVal;
	} // for

	for (unsigned i = 0; i < bSize; i++ ) {
		for (unsigned j = 0; j < pSize; j++) pCond(j, b[i]) -= minVal;
	} // for

	stepNumber = 4;
} // stepSix()
