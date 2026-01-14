#ifndef _APN_
#define _APN_

#include "../common.hpp"
#include "../core/include.hpp"
#include "../statistics/include.hpp"
#include "../ccz/include.hpp"
#include "./ortho_derivative.hpp"


// !SECTION! Invariants themselves

cpp_Spectrum cpp_sigma_multiplicities(
    const cpp_S_box &f,
    const Integer k,
    const Integer n,
    const Integer n_threads);



// !SECTION! Aggregated invariants

std::string cpp_apn_ea_mugshot(
    const cpp_Spectrum &abs_walsh_spec,
    const cpp_Spectrum &deg_spec,
    const cpp_Spectrum &sig_mult,
    const cpp_Spectrum &thk_spec
    );


std::string cpp_apn_ea_mugshot(
    const cpp_S_box & s,
    const unsigned int n_threads
    );
    
#endif
