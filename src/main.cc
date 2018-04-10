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

	unsigned numberOfTrials = 21;
	std::vector<double> variable(numberOfTrials);
	for (unsigned i = 0; i < numberOfTrials; i++) variable[i] = 0.05*(i);
	variable[20] = 0.99;

	unsigned numberOfSimulationsPerTrial = 5000;
	std::vector<std::vector<std::vector<ColVector<double>>>> ospa(numberOfTrials);
	
	for (unsigned i = 0; i < numberOfTrials; i++) {
		std::cout << "\nTrial " << i << std::endl;
		ospa[i].resize(numberOfSimulationsPerTrial);
		for (unsigned j = 0; j < numberOfSimulationsPerTrial; j++) {
			std::cout << "Trial " << i << ", Variable: " << variable[i] << ", Simulation: " << j << std::endl;
			// Declare model;
			rcptr<LinearModel> model = uniqptr<LinearModel>(new LinearModel(20, variable[i]));
			
			// Load the ground truth
			std::vector<std::vector<rcptr<filters::gmm>>> groundTruthBeliefs = model->getGroundTruthBeliefs();

			// Run the filter
			std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates = runLinearGaussianFilter(model);

			// Calculate the OSPA
			ospa[i][j] = calculateOspa(model, groundTruthBeliefs, stateEstimates);

			if (i == 0 && j == 0) {
				// Load ground truth trajectories
				std::vector<std::vector<ColVector<double>>> groundTruthTrajectories = model->getIndividualGroundTruthTrajectories();
				// Get cardinality
				std::vector<unsigned> cardinality = model->getCardinality();
				// Load the measurements
				std::vector<std::vector<ColVector<double>>> measurements = model->getMeasurements();
				// Output results
				outputResults(model, groundTruthTrajectories, measurements, stateEstimates, ospa[i][j], cardinality);
			}
		} // for
	} // for

	// Output results
	printOspaTrials(ospa);
	
	
	//outputResults(model, model->getIndividualGroundTruthTrajectories(), measurements, stateEstimates, ospa, cardinality);
	//std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates = runLinearGaussianFilterMT(model);

	return 0;
}
