/*************************************************************************
 *  Compilation: ./run_main 0
 *  Execution: ./run_main 0
 *  Dependencies:
 *
 * Main app, runs everything.
 *************************************************************************/

// Filter headers
#include <iostream>
#include "system_constants.hpp"
#include "model_declaration.hpp"
#include "algorithmic_steps.hpp"

/**
 * Main app, runs small examples for now.
 *
 * @author SCJ Robertson
 * @since 04/03/18
 */
int main(int, char *argv[]) {
	
	rcptr<LinearModel> model = uniqptr<LinearModel>(new LinearModel());
	std::vector<std::vector<rcptr<filters::gaussian>>> groundTruth = model->getGroundTruthBeliefs();
	std::vector<std::vector<ColVector<double>>> measurements = model->getMeasurements();

	std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates = runLinearGaussianFilter(model);

	return 0;
}
