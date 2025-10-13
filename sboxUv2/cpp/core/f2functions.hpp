#ifndef _f2_funcs_
#define _f2_funcs_

#include "../sboxU.hpp"


// !SECTION! Bit-fiddling


/* The functions defined in this section rely a lot on g++ built-ins
 * intended for bit-fiddling. More info can be found online:
 *
 * https://gcc.gnu.org/onlinedocs/gcc/Bit-Operation-Builtins.html
 *
 */

inline BinWord cpp_msb(const BinWord x)
{
    if (x == 0)
        return 0;
    else
        return 8*(sizeof(BinWord)) - __builtin_clzll(x) - 1;
};

inline BinWord cpp_lsb(const BinWord x)
{
    return __builtin_ctzll(x);
};


inline BinWord cpp_hamming_weight (const BinWord x)
{
    return __builtin_popcountll(x) ;
};


inline BinWord cpp_scal_prod (const BinWord x, const BinWord y)
{
    return __builtin_parity(x & y);
};


inline BinWord cpp_oplus (const BinWord x, const BinWord y)
{
    return x ^ y;
};

BinWord cpp_circ_shift(const BinWord x, int n, int shift);

// !SUBSECTION! tobin and frombin

std::vector<int> cpp_to_bin(const BinWord x, int n);
BinWord cpp_from_bin(const std::vector<int>& v);

// !SECTION! Linear combinations of vectors and their ranks 

std::vector<BinWord> cpp_transpose (const std::vector<BinWord> & l);

BinWord cpp_linear_combination (const std::vector<BinWord> & l, BinWord mask);

Integer cpp_rank_of_vector_set(std::vector<BinWord> l);



// !SECTION! Efficiently enumerating elements in a hyperplane 

/* x: iteration variable name
 * a: difference (BinWord), can be 0, can be greater than size
 * size: must be a power of 2
 * Enumerate in order the x such that all (x,x^a) are distinct
 * Guarantee x < size and x < x^a
 * */
#define FOR_ENUMERATE_DIFFERENCE_COSETS(x,a,size) \
    for(BinWord x = 0, _max = size, _msb = (a == 0 || a > _max) ? _max : 1l << cpp_msb(a); x < _max; x+=_msb)\
        for(BinWord _ceil = x + _msb; x < _ceil; x++)

#endif
