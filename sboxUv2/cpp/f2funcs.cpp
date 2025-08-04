#include "sboxU.hpp"


/* The functions defined in this file rely a lot on g++ built-ins
 * intended for bit-fiddling. More info can be found online:
 *
 * https://gcc.gnu.org/onlinedocs/gcc/Bit-Operation-Builtins.html
 *
 */

uint64_t cpp_msb(uint64_t x)
{
    return 8*(sizeof(BinWord)) - __builtin_clzll(x);
}


uint64_t cpp_lsb(uint64_t x)
{
    return 8*(sizeof(BinWord)) - __builtin_ctzll(x);
}


uint64_t cpp_hamming_weight(uint64_t x)
{
    return __builtin_popcountll(x);
}


uint64_t cpp_scalar_prod (uint64_t x, uint64_t y)
{
    return __builtin_popcountll(x & y);
}


uint64_t cpp_oplus (uint64_t x, uint64_t y)
{
    return x ^ y;
}
