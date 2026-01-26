#ifndef _APN_SN_
#define _APN_SN_


#include "../common.hpp"
#include "../core/s_box.hpp"
#include "../algorithms/include.hpp"
#include "../statistics/differential.hpp"



// Basic libraries
#include <random>




/////////////////
//// TOOLBOX ////
/////////////////

uint64_t switching_neighbors_rand_int_64_cpp();

std::vector<uint64_t> batch_rand_int_cpp(uint64_t size, uint64_t n);

/////////////////////////////
//// SWITCHING NEIGHBORS ////
/////////////////////////////

std::vector<BinWord> cpp_from_vector( Bytearray  v);
BinWord cpp_diff(std::vector<BinWord> f, BinWord a, BinWord x);
std::vector<BinWord> cpp_to_lut_coordinate(std::vector<BinWord> s);
std::vector<BinWord> cpp_complete_basis(std::vector<BinWord> start_basis, int n);
void cpp_sn_add_equations(cpp_S_box f, std::vector<cpp_F2LinearSystem>& E, std::vector<BinWord> indices, uint64_t n_add_eq);
std::vector<std::vector<cpp_S_box>> cpp_non_trivial_sn (cpp_S_box f, uint64_t n_eq, uint64_t n_step);

#endif
