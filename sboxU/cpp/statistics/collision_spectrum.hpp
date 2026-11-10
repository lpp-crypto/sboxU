#ifndef _COLLISION_SPECTRUM_
#define _COLLISION_SPECTRUM_

#include "../common.hpp"
#include "../core/s_box.hpp"
#include "../core/spectrum.hpp"

inline std::vector<Integer> cpp_collision_table(const cpp_S_box &s) {
  std::vector<Integer> result(s.output_space_size(), 0);
#pragma omp parallel for
  for (unsigned int x = 0; x < s.input_space_size(); x++) {
#pragma omp atomic update
    result[s[x]]++;
  }
  return result;
}

inline cpp_Spectrum cpp_f2_collision_spectrum(const cpp_S_box &s) {
  cpp_Spectrum res{};
  const std::vector<Integer> coll_table = cpp_collision_table(s);
  res.incr_by_counting(coll_table);
  return res;
}

#endif
