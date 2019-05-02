/* Time-stamp: <2019-05-02 17:28:33 lperrin>
 *
 * LICENSE
 */ 


#ifndef _SBOXU_CPP_DIFF_LIN_H_
#define _SBOXU_CPP_DIFF_LIN_H_

#include "sboxu_cpp.hpp"
using namespace boost::python;



// !SECTION! Differential

/* @return a python list of python lists d such that
 *
 * d[a][b] = #{x, l[x+a] + l[x] = b},
 *
 * where + denotes XOR.
 * 
 * Throws an error if the length of l is not a power of 2.
 */ 
list ddt(const list& l);


/* @return a python dictionnary s such that
 *
 * s[k] = #{(a, b), ddt[a][b] = k}.
 *
 * The computation is performed using `n_threads` distinct threads.
 *
 * Throws an error if the length of l is not a power of 2.
 */ 
dict differential_spectrum_fast(const list& l, const unsigned int n_threads);


// !SECTION! Linear

Sbox invert_lat_cpp(const std::vector<std::vector<Integer> > l, const unsigned int n);

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


/* @return a python list of python lists w such that
 *
 * w[a][b] = \sum_{x}(-1)^{a.x + b.l[x]},
 *
 * where + denotes XOR, ^ the exponentiation and . the binary scalar
 * product.
 * 
 * Throws an error if the length of l is not a power of 2.
 */
list lat(const list& l);

/* @return a python dictionnary s such that
 *
 * s[k] = #{(a, b), lat[a][b] = k}.
 *
 * The computation is performed using `n_threads` distinct threads.
 *
 * Throws an error if the length of l is not a power of 2.
 */ 
dict walsh_spectrum_fast(const list& l, const unsigned int n_threads);


/** @return a vector containing BinWord (a << n) | b such that
 *
 * L[a, b] == 0,
 *
 * where L is the LAT of s.
 */
std::vector<BinWord> lat_zeroes_cpp(
    const Sbox s,
    const unsigned int n,
    const unsigned int n_threads);


/** @return a vector containing all BinWord:s a such that there exists
 * b where L[a, b] == 0, where L is the LAT of s.
 */
std::vector<BinWord> projected_lat_zeroes_cpp(
    const Sbox s,
    const unsigned int n_threads);



// !SECTION! Boomerang

/* @return a python list of python lists d such that
 *
 * d[a][b] = #{x, S^-1(S(x)+b) + S^-1(S(x+a)+b) = a}
 *
 * where + denotes XOR.
 * 
 * Throws an error if the length of l is not a power of 2.
 */ 
list bct(const list& l);


/* @return a python dictionnary s such that
 *
 * s[k] = #{(a, b), a != 0, b != 0, S^-1(S(x)+b) + S^-1(S(x+a)+b) = a} has k solutions}
 *
 * where + denotes XOR. The computation is performed using `n_threads`
 * distinct threads.
 *
 * Throws an error if the length of l is not a power of 2.
 */ 
dict bct_spectrum_fast(const list& l, const unsigned int n_threads);




#endif /* _SBOXU_CPP_DIFF_LIN_H_ */
