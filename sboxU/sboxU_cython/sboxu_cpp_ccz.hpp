/* Time-stamp: <2020-09-04 15:22:29 lperrin>
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



#endif /* _SBOXU_CPP_CCZ_H_ */