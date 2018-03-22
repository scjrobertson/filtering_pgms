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
	
	return 0.0;
} // hungarianCost()

void hungarianStepOne(Matrix<double> & pCond,
		unsigned *stepNumber) {
	unsigned numberOfRows = pCond.rows();
	unsigned numberOfColumns = pCond.cols();

	for (unsigned i = 0; i < numberOfRows; i++) {
		double minValue = std::numeric_limits<double>::infinity();
		for (unsigned j = 0; j < numberOfColumns; j++) minValue = std::min(minValue, pCond(i, j));
		for (unsigned j = 0; j < numberOfColumns; j++) pCond(i, j) = pCond(i,j) - minValue;
	} // for

	*stepNumber = 1;
} // stepOne()

void hungarianStepTwo(Matrix<double> & pCond,
		Matrix<unsigned> & mask,
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		unsigned *stepNumber) {
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

	*stepNumber = 3;
} // stepTwo()

void hungarianStepThree(Matrix<unsigned> & mask,
		ColVector<unsigned> & columnsCovered,
		unsigned *stepNumber) {
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

	if (totalSum == numberOfColumns) *stepNumber = 7;
	else *stepNumber = 4;
} // stepThree()

void hungarianStepFour(Matrix<double> & pCond, 
		Matrix<unsigned> & mask, 
		ColVector<unsigned> & rowsCovered,
		ColVector<unsigned> & columnsCovered,
		unsigned *stepNumber) {

	unsigned pSize = pCond.cols();
	std::vector<unsigned> zRow, zCol;
	bool zFlag = true;

	while(zFlag) {
		unsigned row = 0; unsigned col = 0;
		bool exitFlag = true;
		unsigned i = 1; unsigned j = 1;

		while (exitFlag) {
			if (pCond(i, j) == 0 && rowsCovered[i] == 0 && columnsCovered[j] == 0 ) {
				row = i; 
				col = j;
				exitFlag = false;
			} // if
			j++;
			if (j > pSize) j = 1; i++;
			if (i > pSize) exitFlag = false;
		} // while

		if (row == 0) {
			*stepNumber = 6;
			zFlag = false;
			zRow.push_back(0); zCol.push_back(0);
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
				for (unsigned i = 0; i < indices.size(); i++) {
					columnsCovered[indices[i]] = 0;
					zCol.push_back(indices[i]);
				} // for

			} else {
				*stepNumber = 5;
				zFlag = false;
				zRow.push_back(row);
				zCol.push_back(col);
			} // if
		} // if
	} // while
} // stepFour()

Matrix<unsigned> stepFive(Matrix<unsigned> mask,
		std::vector<unsigned> zRows,
		std::vector<unsigned> zCols,
		ColVector<unsigned> rowsCovered,
		ColVector<unsigned> columnsCovered,
		unsigned *stepNumber) {
	
	unsigned pSize = mask.cols();

	bool zFlag = true;
	unsigned i = 1;
	while (zFlag) {
		unsigned rowIndex = -1;
		for (unsigned j = 0; j < pSize; j++) if ( mask(j, zCols[i]) == 1 ) rowIndex = j;

		if (rowIndex > 0) {
			i++;
			zRows[i] = rowIndex;
			zCols[i] = zCols[i-1];
		} else {
			zFlag = false;
		} // if
		
		if (zFlag == true) {
			unsigned columnIndex = -1;
			for (unsigned j = 0; j < pSize; j++) if ( mask(zRows[i], j) == 2 ) columnIndex = j;
			i++;
			zRows[i] = zRows[i-1];
			zCols[i] = columnIndex;
		} // if
	} // while

	for (unsigned i = 0; i < zRows.size(); i++) {
		if (mask(zRows[i], zCols[i]) == 1 ) mask(zRows[i], zCols[i]) = 1;
		else mask(zRows[i], zCols[i]) = 0;
	} // for

	for (unsigned i = 0; i < zRows.size(); i++) zRows[i] = 0;
	for (unsigned i = 0; i < zCols.size(); i++) zCols[i] = 0;

	for (unsigned i = 0; i < pSize; i++) {
		for (unsigned j = 0; j < pSize; j++) if (mask(i, j) == 2) mask(i, j) = 0;
	}

	return mask;
} // stepFive()

Matrix<double> stepSix(Matrix<double> pCond,
		ColVector<unsigned> rowsCovered,
		ColVector<unsigned> columnsCovered,
		unsigned *stepNumber) {
	
	std::vector<unsigned> a, b, c, d;
	for (unsigned i = 0; i < rowsCovered.getNElem(); i++) {
		if(rowsCovered[i] == 0) a.push_back(i);
		if(rowsCovered[i] == 1) d.push_back(i);
	} // for
	for (unsigned i = 0; i < columnsCovered.getNElem(); i++) {
		if(columnsCovered[i] == 0) b.push_back(i);
		if(columnsCovered[i] == 1) d.push_back(i);
	} // for

	double minVal = std::numeric_limits<double>::infinity();

	unsigned rSize = rowsCovered.size();
	unsigned cSize = columnsCovered.size();

	for (unsigned i = 0; i < rSize; i++) {
		for (unsigned j = 0; j < cSize; j++) {
			minVal = std::min(minVal, pCond(a[i], b[j]));
		} // for
	} // for

	unsigned pSize = pCond.cols();
	unsigned bSize = b.size(); 
	unsigned dSize = d.size();
	
	for (unsigned i = 0; i < dSize; i++ ) {
		for (unsigned j = 0; j < pSize; j++) pCond(d[i], j) += minVal;
	} // for

	for (unsigned i = 0; i < bSize; i++ ) {
		for (unsigned j = 0; j < pSize; j++) pCond(j, b[i]) -= minVal;
	} // for

	return pCond;
} // stepSix()
