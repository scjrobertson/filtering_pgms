#ifndef ALGORITHMIC_STEPS_HPP
#define ALGORITHMIC_STEPS_HPP

#include <iostream>
#include "system_constants.hpp"
#include "model_declaration.hpp"
#include "utils.hpp"

// Standard Filter
std::vector<std::vector<rcptr<filters::gmm>>> runLinearGaussianFilter(rcptr<LinearModel> model);

std::vector<rcptr<filters::gmm>> predictMultipleTargetsLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> targets);

std::vector<rcptr<filters::updateComponents>> createMultipleUpdateComponentsLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> targets);

std::vector<std::vector<rcptr<filters::gmm>>> createUpdateOptionsLinear(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> predictedStates,
		std::vector<rcptr<filters::updateComponents>> kalmanComponents,
		std::vector<ColVector<double>> measurements);

Matrix<double> createAssociationMatrix(rcptr<LinearModel> model,
		unsigned numberOfMeasurements,
		std::vector<std::vector<rcptr<filters::gmm>>> updateOptions);

Matrix<double> loopyBeliefUpdatePropagation(rcptr<LinearModel> model,
		Matrix<double> associationMatrix);

std::vector<rcptr<filters::gmm>> updateTargetStatesLinear(rcptr<LinearModel> model,
		std::vector<std::vector<rcptr<filters::gmm>>> updateOptions,
		Matrix<double> associationMatrix,
		Matrix<double> updatedAssociations);

// Measurement-to-Target Filter
std::vector<std::vector<rcptr<filters::gmm>>> runLinearGaussianFilterMT(rcptr<LinearModel> model);

std::vector<rcptr<filters::cfm>> createCanonicalLikelihoods(rcptr<LinearModel> model,
		std::vector<ColVector<double>> measurements);

std::vector<rcptr<filters::gmm>> updateTargetStatesLinearMT(rcptr<LinearModel> model,
		std::vector<rcptr<filters::gmm>> predictedStates,
		std::vector<rcptr<filters::cfm>> likelihoods,
		Matrix<double> updatedAssociations);

std::vector<std::vector<rcptr<filters::gmm>>> createUpdateOptionsLinearMT(rcptr<LinearModel> model, 
		std::vector<rcptr<filters::gmm>> predictedStates,
		std::vector<rcptr<filters::updateComponents>> kalmanComponents,
		std::vector<ColVector<double>> measurements);

#endif // ALGORITHMIC_STEPS_HPP
