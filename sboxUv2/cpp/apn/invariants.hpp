#ifndef _APN_
#define _APN_

#include "../sboxU.hpp"
#include "../core/include.hpp"
#include "../statistics/include.hpp"


cpp_S_box cpp_ortho_derivative(const cpp_S_box q);

cpp_Spectrum cpp_sigma_multiplicities(
    const cpp_S_box f,
    const Integer k,
    const Integer n,
    const Integer n_threads);

    
#endif
