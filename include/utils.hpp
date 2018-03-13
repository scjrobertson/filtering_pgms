#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include "emdw.hpp"
#include "genvec.hpp"
#include "genmat.hpp"
#include "system_constants.hpp"

rcptr<filters::gmm> weakMarginalisation (rcptr<filters::gmm> gmm);

#endif // UTILS_HPP
