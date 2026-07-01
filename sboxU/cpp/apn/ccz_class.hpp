#ifndef _APN_CCZ_CLASS_
#define _APN_CCZ_CLASS_


#include "../common.hpp"
#include "../core/include.hpp"
#include "../ccz/include.hpp"
#include "./invariants.hpp"
#include "./ortho_derivative.hpp"


std::vector<cpp_F2AffineMap> cpp_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads,
    const std::string & mode = "standard"
    );

std::vector<cpp_F2AffineMap> cpp_ea_mappings_from_ortho_derivative(
    const cpp_S_box & s,
    const cpp_S_box & s_prime,
    const unsigned int n_threads
    );
    
std::vector<cpp_S_box> cpp_enumerate_ea_classes_quadratic_apn(
    const cpp_S_box &s,
    const unsigned int n_threads,
    const std::string & mode = "standard"
    );


cpp_S_box cpp_ccz_equivalent_quadratic_function(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

std::vector<cpp_F2AffineMap>cpp_graph_automorphisms_from_derivatives(cpp_S_box s);

std::vector<cpp_F2AffineMap> cpp_gen_set_graph_automorphisms_from_derivative(
    const cpp_S_box & s
    );


std::vector<cpp_F2AffineMap> cpp_graph_el_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

std::vector<cpp_F2AffineMap> cpp_gen_set_F2AffineMap_group(
    std::vector<cpp_F2AffineMap> G,
    const std::string & mode = "deterministic"
    );

#endif
