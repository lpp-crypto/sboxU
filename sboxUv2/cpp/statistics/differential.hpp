#ifndef _STAT_DIFFERENTIAL_
#define _STAT_DIFFERENTIAL_

#include "../sboxU.hpp"
#include "../s_box.hpp"


cpp_Spectrum cpp_differential_spectrum(
    const cpp_S_box s,
    const unsigned int n_threads
    );

std::vector<Integer> cpp_ddt_row(const cpp_S_box s, const BinWord delta);

std::vector< std::vector<Integer> > cpp_ddt(const cpp_S_box s) ;

bool cpp_is_differential_uniformity_smaller_than(
    const cpp_S_box s,
    const Integer u
    ) ;


#endif
