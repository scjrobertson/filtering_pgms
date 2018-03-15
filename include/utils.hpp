#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include "emdw.hpp"
#include "genvec.hpp"
#include "genmat.hpp"
#include "system_constants.hpp"

rcptr<filters::gmm> weakMarginalisation (rcptr<filters::gmm> gmm);

bool haveIntersectingDomains(std::vector<unsigned short> a, 
		std::vector<unsigned short> b);

#endif // UTILS_HPP
