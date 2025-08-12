#ifndef _ALGO_SP_SEARCH
#define _ALGO_SP_SEARCH

#include "../sboxU.hpp"


std::vector<BinWord> cpp_extract_vector(
    const std::vector<BinWord> & z,
    const BinWord a
    );


std::vector<std::vector<BinWord> > cpp_extract_bases(
    std::vector<BinWord> & z,
    const Integer dimension,
    Integer n_threads,
    const std::string end_condition
    );


std::vector<std::vector<BinWord> > cpp_extract_affine_bases(
    std::vector<BinWord> & z,
    const Integer dimension,
    Integer n_threads,
    const std::string end_condition
    );



#endif
