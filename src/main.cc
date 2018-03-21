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
#include "ospa.hpp"

/**
 * Main app, runs small examples for now.
 *
 * @author SCJ Robertson
 * @since 04/03/18
 */
int main(int, char *argv[]) {
	
	rcptr<LinearModel> model = uniqptr<LinearModel>(new LinearModel());
	std::vector<std::vector<rcptr<filters::gmm>>> groundTruthBeliefs = model->getGroundTruthBeliefs();
	std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates = runLinearGaussianFilter(model);

	ColVector<double> ospa = calculateOspa(groundTruthBeliefs[15], stateEstimates[4], model->ospaC, model->ospaP);

	return 0;
}
