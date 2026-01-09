#ifndef _ANF_HPP_
#define _ANF_HPP_

#include "s_box.hpp"
#include "../sboxU.hpp"
#include "f2functions.hpp"
#include "spectrum.hpp"


std::vector<BinWord> cpp_anf_component( const cpp_S_box &f);

std::vector<BinWord>cpp_quadratic_compact_representation( const cpp_S_box &f);

std::vector<BinWord> cpp_quadratic_sbox_from_compact_representation( std::vector<BinWord> compact_representation, int64_t n, int64_t m);

Integer cpp_degree_component(const cpp_S_box &f);

cpp_Spectrum cpp_degree_spectrum(const cpp_S_box &f);

cpp_Spectrum cpp_monomial_degree_spectrum_component(const cpp_S_box &f);

Integer cpp_algebraic_degree(const cpp_S_box &f);
    
#endif
