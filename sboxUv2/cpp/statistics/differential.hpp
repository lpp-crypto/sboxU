#ifndef _STAT_DIFFERENTIAL_
#define _STAT_DIFFERENTIAL_

#include "../common.hpp"
#include "../core/include.hpp"


cpp_Spectrum cpp_differential_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

cpp_Spectrum cpp_differential_spectrum_fast(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

std::vector<Integer> cpp_ddt_row(
    const cpp_S_box & s,
    const BinWord delta
    );

std::vector< std::vector<Integer> > cpp_ddt(const cpp_S_box & s) ;


std::vector< std::vector<BinWord>> cpp_xddt_row(const cpp_S_box & s, const BinWord delta);

Xtable cpp_xddt(const cpp_S_box & s);

std::vector<BinWord> cpp_xddt_entry(const cpp_S_box &s, const BinWord a, const BinWord b);


std::vector< std::vector<BinWord>> cpp_yddt_row(const cpp_S_box & s, const BinWord delta);
Xtable cpp_yddt(const cpp_S_box & s);


std::vector< std::vector<BinWord>> cpp_zddt_row(const cpp_S_box & s, const BinWord delta);
Xtable cpp_zddt(const cpp_S_box & s);

bool cpp_is_differential_uniformity_smaller_than(
    const cpp_S_box & s,
    const Integer u
    ) ;


#endif
