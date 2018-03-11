#ifndef ALGORITHMIC_STEPS_HPP
#define ALGORITHMIC_STEPS_HPP

#include "system_constants.hpp"
#include "model_declaration.hpp"

std::vector<std::vector<rcptr<filters::gmm>>> runLinearGaussianFilter(rcptr<LinearModel> model);



#endif // ALGORITHMIC_STEPS_HPP
