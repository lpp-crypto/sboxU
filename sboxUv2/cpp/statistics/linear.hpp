#ifndef _STAT_LINEAR_
#define _STAT_LINEAR_


#include "../sboxU.hpp"
#include "../s_box.hpp"


std::vector<Integer> cpp_walsh_transform(const cpp_S_box f);

cpp_Spectrum cpp_walsh_spectrum(cpp_S_box s);

std::vector< std::vector<Integer> > cpp_lat(cpp_S_box s);

cpp_S_box cpp_invert_lat(const std::vector< std::vector<Integer> > l);


#endif
