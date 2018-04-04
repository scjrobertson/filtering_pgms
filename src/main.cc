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

	unsigned numberOfTrials = 10;
	std::vector<unsigned> variable(numberOfTrials);
	for (unsigned i = 0; i < numberOfTrials; i++) variable[i] = 5*(i+1);

	unsigned numberOfSimulationsPerTrial = 1000;
	std::vector<std::vector<std::vector<ColVector<double>>>> ospa(numberOfTrials);
	
	for (unsigned i = 0; i < numberOfTrials; i++) {
		ospa[i].resize(numberOfSimulationsPerTrial);
		for (unsigned j = 0; j < numberOfSimulationsPerTrial; j++) {
			// Declare model;
			rcptr<LinearModel> model = uniqptr<LinearModel>(new LinearModel(variable[i]));
			
			// Load the ground truth
			std::vector<std::vector<rcptr<filters::gmm>>> groundTruthBeliefs = model->getGroundTruthBeliefs();
			std::vector<unsigned> cardinality = model->getCardinality();

			// Load the measurements
			std::vector<std::vector<ColVector<double>>> measurements = model->getMeasurements();
			
			// Run the filter
			std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates = runLinearGaussianFilter(model);

			// Calculate the OSPA
			ospa[i][j] = calculateOspa(model, groundTruthBeliefs, stateEstimates);
		} // for
	} // for

	// Output results
	printOspaTrials(ospa);
	
	
	//outputResults(model, model->getIndividualGroundTruthTrajectories(), measurements, stateEstimates, ospa, cardinality);
	//std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates = runLinearGaussianFilterMT(model);

	return 0;
}
