#ifndef _f2_funcs_
#define _f2_funcs_

#include "sboxU.hpp"


// !SECTION! Bit-fiddling

BinWord cpp_msb (BinWord x);

BinWord cpp_lsb (BinWord x);

BinWord cpp_hamming_weight (BinWord x);

BinWord cpp_scal_prod (BinWord x, BinWord y);

BinWord cpp_oplus (BinWord x, BinWord y);


// !SECTION! Linear combinations of vectors and their ranks 

BinWord cpp_linear_combination (std::vector<BinWord> l, BinWord mask);

Integer cpp_rank_of_vector_set(std::vector<BinWord> l);


class cpp_Linear_basis
{
private:
    std::map<Integer, BinWord> basis;
public:
    cpp_Linear_basis(std::vector<BinWord> l) ;

    void add_to_span(BinWord x);

    std::vector<BinWord> get_basis() const;

    Integer rank() const;
    
    std::vector<BinWord> span() const;
};
     
#endif
