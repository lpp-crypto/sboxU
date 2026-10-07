#ifndef _INVARIANTS_
#define _INVARIANTS_

#include "../common.hpp"
#include "../core/include.hpp"
#include "./bigvectors.hpp"

std::vector<cpp_BigF2Vector> cpp_basis_invariants_from_cycles(const cpp_S_box S);
std::vector<BinWord> cpp_vectors_of_hamming_weight(BinWord h, BinWord n);
std::tuple<Lut, BinWord> yann_permutation(BinWord d, BinWord n);
std::vector<cpp_BigF2Vector> cpp_all_invariants_up_to_degree(const cpp_S_box S, BinWord d);

#endif