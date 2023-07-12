#ifndef _SBOXU_CPP_FP_LAT_H_
#define _SBOXU_CPP_FP_LAT_H_

#include "sboxu_cpp.hpp"
#include <pocketfft_hdronly.h>
#include <omp.h>

std::vector<std::vector<double>> fpt_lat(const std::vector<int>& s, const int p, const int m, const unsigned int num_threads);
std::vector<double> fpt_lat_row(const std::vector<int>& s, const int p, const int m, const int a);
std::vector<double> fpt_lat_column(const std::vector<int>& s, const int p, const int m, const int b);
double fpt_max_lat(const std::vector<int>& s, const int p, const int m, const unsigned int num_threads);
std::map<double, int> fpt_walsh_spectrum(const std::vector<int>& s, const int p, const int m, const double epsilon, const unsigned int num_threads);

#endif
