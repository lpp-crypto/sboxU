#ifndef _APN_CCZ_CLASS_
#define _APN_CCZ_CLASS_


#include "../common.hpp"
#include "../core/include.hpp"
#include "../ccz/include.hpp"
#include "./invariants.hpp"


std::vector<cpp_BinLinearMap> cpp_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

std::vector<cpp_BinLinearMap> cpp_ea_mappings_from_ortho_derivative(
    const cpp_S_box & s,
    const cpp_S_box & s_prime,
    const unsigned int n_threads
    );
    
std::vector<cpp_S_box> cpp_enumerate_ea_classes_quadratic_apn(
    const cpp_S_box &s,
    const unsigned int n_threads
    );


cpp_S_box cpp_ccz_equivalent_quadratic_function(
    const cpp_S_box & s,
    const unsigned int n_threads
    );


#endif
