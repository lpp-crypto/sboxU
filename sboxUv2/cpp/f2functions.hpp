#ifndef _f2_funcs_
#define _f2_funcs_

#include "sboxU.hpp"

BinWord cpp_msb (BinWord x);

BinWord cpp_lsb (BinWord x);

BinWord cpp_hamming_weight (BinWord x);

BinWord cpp_scalar_prod (BinWord x, BinWord y);

BinWord cpp_oplus (BinWord x, BinWord y);

BinWord cpp_linear_combination (std::vector<BinWord> v, BinWord mask);
     
#endif
