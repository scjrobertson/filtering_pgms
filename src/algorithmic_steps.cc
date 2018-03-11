#include "algorithmic_steps.hpp"

std::vector<std::vector<rcptr<filters::gmm>>> runLinearGaussianFilter(rcptr<LinearModel> model) {
	std::vector<std::vector<rcptr<filters::gmm>>> targets; targets.clear(); targets.resize(model->simulationLength);

	for (unsigned i = 0; i < model->simulationLength; i++) {
		std::cout << i << std::endl;
	} // for

	return targets;
} // runLinearGaussianFilter()
