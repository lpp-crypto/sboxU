/* Time-stamp: <2021-09-21 11:09:03 lperrin>
 *
 * LICENSE
 */ 


#ifndef _SBOXU_CPP_DIFF_LIN_H_
#define _SBOXU_CPP_DIFF_LIN_H_

#include "sboxu_cpp.hpp"
#include "sboxu_cpp_utils.hpp"


// !SECTION! Differential

/*  @return a c++ vector of c++ vectors v such that
 *
 * v[a][b] = #{x, s[x+a] + s[x] = b},
 *
 * where + denotes XOR.
 * 
 * Throws an error if the length of s is not a power of 2. Not yet
 */ 
std::vector< std::vector<Integer> > ddt_cpp(const Sbox s);

/* @return a c++ map m such that
 *
 * m[k] = #{(a, b), ddt[a][b] = k}.
 *
 * The computation is performed using `n_threads` distinct threads.
 *
 * Throws an error if the length of s is not a power of 2. Not yet
 */ 
std::map<Integer,Integer> differential_spectrum_fast(const Sbox s, const unsigned int n_threads);


/* @return true if and only if the differential uniformity of s is at
 * most equal to u.
 */ 
bool is_differential_uniformity_smaller_than_cpp(const Sbox  s, const Integer u);


/**@return a vector v such that v[c] is a map of the form {k : n_k}
 * such that s[x^a] ^ c s[x] = b has exactly k solutions x for n_k
 * different pairs (a,b).
 */
std::vector<std::map<Integer, Integer> > c_differential_spectra_hpp(
    const Sbox s,
    const Sbox l_table,
    const Sbox e_table);

// !SECTION! Linear

Sbox invert_lat_cpp(const std::vector< std::vector<Integer> > l, const unsigned int n);

/* @return the Fourier transform of the Boolean function f, i.e. an
 * array T such that
 *
 * T(a) = \\sum_{x \\in \\{ 0,1 \\}^n}} (-1)^{<a,x> + f(x)} .
 *
 * f is assumed to be the truth-table of a Boolean function. Its size
 * must thus be a power of two. At the end of the function, the output
 * T is the Walsh spectrum of f. If n is the size of the input of f,
 * then:
 */
std::vector<Integer> walsh_spectrum_coord(const Sbox f);

/*  @return a c++ vector of c++ vectors v such that
 *
 * w[a][b] = \sum_{x}(-1)^{a.x + b.l[x]},
 *
 * where + denotes XOR, ^ the exponentiation and . the binary scalar
 * product.
 * 
 * Throws an error if the length of l is not a power of 2.
 */
std::vector< std::vector<Integer> > lat_cpp(const Sbox s);

/* @return a c++ map count such that
 *
 * count[k] = #{(a, b), lat[a][b] = k}.
 *
 * The computation is performed using `n_threads` distinct threads.
 *
 * Throws an error if the length of l is not a power of 2.
 */ 
std::map<Integer, Integer> walsh_spectrum_fast_cpp(const Sbox s, const unsigned int n_threads);

/** @return a vector containing BinWord (a << n) | b such that
 *
 * L[a, b] == 0,
 *
 * where L is the LAT of s.
 */
std::vector<BinWord> lat_zeroes_cpp(const Sbox s, const unsigned int n, const unsigned int n_threads);


/** @return a vector containing all BinWord:s a such that there exists
 * b where L[a, b] == 0, where L is the LAT of s.
 */
std::vector<BinWord> projected_lat_zeroes_cpp( const Sbox s, const unsigned int n_threads);


// !SECTION! Boomerang

/* @return a c++ vector of c++ vectors v such that
 *
 * v[a][b] = #{x, S^-1(S(x)+b) + S^-1(S(x+a)+b) = a}
 *
 * where + denotes XOR.
 * 
 * Throws an error if the length of l is not a power of 2.
 */ 
std::vector< std::vector<Integer>> bct_cpp(const Sbox s);

/* @return a c++ map s such that
 *
 * s[k] = #{(a, b), a != 0, b != 0, S^-1(S(x)+b) + S^-1(S(x+a)+b) = a} has k solutions}
 *
 * where + denotes XOR. The computation is performed using `n_threads`
 * distinct threads.
 *
 * Throws an error if the length of l is not a power of 2.
 */ 
std::map<Integer, Integer> bct_spectrum_fast_cpp(const Sbox s, const unsigned int n_threads);


// !SECTION! Quadratic functions

/* @return an Sbox containing the LUT of the ortho-derivative of
 * the function s.
 *
 * If s is not crooked or if it is not APN, returns an empty Sbox.
 */ 
Sbox ortho_derivative_fast(const Sbox& s);

#endif /* _SBOXU_CPP_DIFF_LIN_H_ */
