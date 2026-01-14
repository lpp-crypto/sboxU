#ifndef _ANF_HPP_
#define _ANF_HPP_

#include "s_box.hpp"
#include "../common.hpp"
#include "f2functions.hpp"
#include "spectrum.hpp"


std::vector<BinWord> cpp_anf_component( const cpp_S_box &f);

Integer cpp_degree_component(const cpp_S_box &f);

cpp_Spectrum cpp_degree_spectrum(const cpp_S_box &f);

cpp_Spectrum cpp_monomial_degree_spectrum_component(const cpp_S_box &f);

Integer cpp_algebraic_degree(const cpp_S_box &f);
    
#endif
