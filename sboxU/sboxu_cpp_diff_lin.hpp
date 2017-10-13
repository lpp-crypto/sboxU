/* Time-stamp: <2017-10-06 19:10:15 lperrin>
 *
 * LICENSE
 */ 


#ifndef _SBOXU_CPP_DIFF_LIN_H_
#define _SBOXU_CPP_DIFF_LIN_H_

#include "sboxu_cpp.hpp"
using namespace boost::python;


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


#endif /* _SBOXU_CPP_DIFF_LIN_H_ */
