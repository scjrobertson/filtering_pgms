#ifndef ALGORITHMIC_STEPS_HPP
#define ALGORITHMIC_STEPS_HPP

#include <iostream>
#include "system_constants.hpp"
#include "model_declaration.hpp"
#include "utils.hpp"

std::vector<std::vector<rcptr<filters::gmm>>> runLinearGaussianFilter(rcptr<LinearModel> model);

std::vector<rcptr<filters::gmm>> predictMultipleTargetsLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> targets);

std::vector<rcptr<filters::updateComponents>> createMultipleUpdateComponentsLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> targets);

std::vector<std::vector<rcptr<filters::gmm>>> updateMultipleLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> predictedStates,
		std::vector<rcptr<filters::updateComponents>> kalmanComponents,
		std::vector<ColVector<double>> measurements);

Matrix<double> createAssociationMatrix(rcptr<LinearModel> model,
		unsigned numberOfMeasurements,
		std::vector<std::vector<rcptr<filters::gmm>>> updateOptions);

Matrix<double> loopyBeliefUpdatePropagation(rcptr<LinearModel> model,
		Matrix<double> associationMatrix);

#endif // ALGORITHMIC_STEPS_HPP
