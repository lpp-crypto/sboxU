#ifndef _STAT_DIFFERENTIAL_
#define _STAT_DIFFERENTIAL_

#include "../sboxU.hpp"
#include "../s_box.hpp"

cpp_Spectrum cpp_differential_spectrum(
    const cpp_S_box s,
    const unsigned int n_threads
    );

std::vector< std::vector<Integer> > cpp_ddt(const cpp_S_box s) ;


#endif
