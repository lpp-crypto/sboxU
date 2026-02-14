#ifndef _QUASI_DIFFERENTIAL_
#define _QUASI_DIFFERENTIAL_

#include "differential.hpp"
#include "../core/include.hpp"

Integer cpp_qddt_coeff(const cpp_S_box &s, const BinWord a, const BinWord b, const BinWord u, const BinWord v);

std::vector<std::vector<Integer>> cpp_qddt_fixed_differential(const cpp_S_box &s, const BinWord a, const BinWord b);


#endif