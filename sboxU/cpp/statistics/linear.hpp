#ifndef _STAT_LINEAR_
#define _STAT_LINEAR_


#include "../common.hpp"
#include "../core/include.hpp"


std::vector<Integer> cpp_walsh_transform(const cpp_S_box & f);

cpp_Spectrum cpp_walsh_spectrum(const cpp_S_box & s, unsigned int n_threads);

cpp_Spectrum cpp_absolute_walsh_spectrum(const cpp_S_box & s, unsigned int n_threads);

std::vector< std::vector<Integer> > cpp_lat(const cpp_S_box & s);

cpp_S_box cpp_invert_lat(const std::vector< std::vector<Integer> > & l);


#endif
