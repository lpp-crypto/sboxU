/* Time-stamp: <2021-08-12 17:10:24 lperrin>
 *
 * LICENSE
 */ 

#ifndef _SBOXU_CPP_CCZ_H_
#define _SBOXU_CPP_CCZ_H_

#include "sboxu_cpp.hpp"



// !SECTION! Python-facing functions

std::vector<BinWord> extract_vector_cpp(
    const std::vector<BinWord> & z,
    const BinWord a);

    
std::vector<std::vector<BinWord> > extract_bases_cpp(
    std::vector<BinWord> & z,
    const Integer dimension,
    const Integer word_length,
    Integer n_threads,
    const std::string end_condition);

std::vector<std::vector<BinWord> > extract_affine_bases_cpp(
    std::vector<BinWord> & z,
    const Integer dimension,
    const Integer word_length,
    Integer n_threads,
    const std::string end_condition);


// !SECTION! Sigma multiplicities

std::map<BinWord, Integer> sigma_multiplicities_cpp(
    const Sbox f,
    const Integer k,
    const Integer n,
    const Integer n_threads);


#endif /* _SBOXU_CPP_CCZ_H_ */
