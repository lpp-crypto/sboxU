#include "sboxU.hpp"


/* The functions defined in this file rely a lot on g++ built-ins
 * intended for bit-fiddling. More info can be found online:
 *
 * https://gcc.gnu.org/onlinedocs/gcc/Bit-Operation-Builtins.html
 *
 */

BinWord cpp_msb(BinWord x)
{
    if (x == 0)
        return 0;
    else
        return 8*(sizeof(BinWord)) - __builtin_clzll(x) - 1;
}


BinWord cpp_lsb(BinWord x)
{
    return __builtin_ctzll(x);
}


BinWord cpp_hamming_weight (BinWord x)
{
    return __builtin_popcountll(x) ;
}


BinWord cpp_scal_prod (BinWord x, BinWord y)
{
    return __builtin_popcountll(x & y);
}


BinWord cpp_oplus (BinWord x, BinWord y)
{
    return x ^ y;
}


BinWord cpp_linear_combination (std::vector<BinWord> v, BinWord mask)
{
    BinWord result = 0;
    for(auto &x : v)
    {
        if (mask & 1)
            result ^= x;
        mask >>= 1;
    }
    return result;
}
