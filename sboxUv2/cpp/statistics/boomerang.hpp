#ifndef _STAT_BOOMERANG_
#define _STAT_BOOMERANG_

#include "../sboxU.hpp"
#include "../core/include.hpp"

cpp_Spectrum cpp_bct_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

std::vector< std::vector<Integer> > cpp_bct(const cpp_S_box & s) ;


#endif
