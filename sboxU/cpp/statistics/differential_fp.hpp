#ifndef _STAT_DIFFERENTIAL_FP_
#define _STAT_DIFFERENTIAL_FP_

#include "../common.hpp"
#include "../core/s_box_fp.hpp"
#include "../core/spectrum.hpp"

/** Computes one row of the DDT for an F_p S-box.
 *
 *  @param s          The S-box.
 *  @param delta_int  The integer index of the input difference.
 *  @return A vector `row` of size p^u such that `row[gamma_int]` equals
 *          #{x : S(x+delta) - S(x) = gamma (mod p)}.
 */
inline std::vector<Integer> cpp_fp_ddt_row(const cpp_S_box_fp &s,
                                           const int delta_int) {
  const int p = s.get_p();
  const std::vector<FpWord> &input_space = s.get_input_space();
  const std::vector<FpWord> &lut = s.get_lut();
  const std::vector<Integer> &powers_in = s.get_powers_in();
  const std::vector<Integer> &powers_out = s.get_powers_out();
  const int n = input_space.size();
  const int m = s.get_output_space().size();
  const int in_dim = s.get_input_size();
  const int out_dim = s.get_output_size();
  const FpWord delta = cpp_S_box_fp::int_to_vec(delta_int, input_space);

  std::vector<Integer> result(m, 0);
  for (int i = 0; i < n; i++) {
    const FpWord &x = input_space[i];
    FpWord x_plus_delta(in_dim);
    for (int j = 0; j < in_dim; j++)
      x_plus_delta[j] = (x[j] + delta[j]) % p;
    Integer xd_int = cpp_S_box_fp::vec_to_int(x_plus_delta, powers_in);
    FpWord gamma(out_dim);
    for (int j = 0; j < out_dim; j++)
      gamma[j] = (lut[xd_int][j] + p - lut[i][j]) % p;
    result[cpp_S_box_fp::vec_to_int(gamma, powers_out)]++;
  }
  return result;
}

/** Computes the full DDT of an F_p S-box.
 *
 *  @param s  The S-box.
 *  @return A 2D array `table` such that `table[delta_int][gamma_int]` equals
 *          #{x : S(x+delta) - S(x) = gamma (mod p)}.
 *          The row for delta=0 has `table[0][0] = p^t` and all others 0.
 */
inline std::vector<std::vector<Integer>> cpp_fp_ddt(const cpp_S_box_fp &s) {
  const int n = s.get_input_space().size();
  const int m = s.get_output_space().size();
  std::vector<std::vector<Integer>> table(n, std::vector<Integer>(m, 0));
  table[0][0] = n;
  for (int delta_int = 1; delta_int < n; delta_int++)
    table[delta_int] = cpp_fp_ddt_row(s, delta_int);
  return table;
}

/** Computes the differential spectrum of an F_p S-box.
 *
 *  Iterates over all non-zero input differences and tallies the DDT
 *  coefficients using OpenMP parallelism. The delta=0 row is excluded
 *  because it is trivially [p^t, 0, ..., 0].
 *
 *  @param s  The S-box.
 *  @return A `cpp_Spectrum` `sp` such that `sp[k]` is the number of pairs
 *          (delta, gamma) with delta != 0 and D[delta][gamma] = k.
 */
inline cpp_Spectrum cpp_fp_differential_spectrum(const cpp_S_box_fp &s) {
  cpp_Spectrum count;
  const int n = s.get_input_space().size();
#pragma omp parallel for reduction(aggregateSpectrum : count)
  for (int delta_int = 1; delta_int < n; delta_int++)
    count.incr_by_counting(cpp_fp_ddt_row(s, delta_int));
  return count;
}

#endif
