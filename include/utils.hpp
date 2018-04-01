#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include "emdw.hpp"
#include "genvec.hpp"
#include "genmat.hpp"
#include "system_constants.hpp"
#include "model_declaration.hpp"

rcptr<filters::gmm> weakMarginalisation(rcptr<filters::gmm> gmm);

rcptr<filters::gmm> gaussianMixturePruning(rcptr<filters::gmm> gmm,
		double componentWeightThreshold,
		double componentUnionDistance,
		unsigned maximumNumberOfMixtureComponents);

bool haveIntersectingDomains(std::vector<unsigned short> a, 
		std::vector<unsigned short> b);

void outputResults(rcptr<LinearModel> model,
		std::vector<std::vector<ColVector<double>>> groundTruth,
		std::vector<std::vector<ColVector<double>>> measurements,
		std::vector<std::vector<rcptr<filters::gmm>>> stateEstimates,
		std::vector<ColVector<double>> ospa,
		std::vector<unsigned> cardinality);

#endif // UTILS_HPP
