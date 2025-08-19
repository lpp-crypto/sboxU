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



// !SECTION! Linear combinations of vectors and their ranks 

std::vector<BinWord> cpp_transpose (const std::vector<BinWord> & l);

BinWord cpp_linear_combination (const std::vector<BinWord> & l, BinWord mask);

Integer cpp_rank_of_vector_set(std::vector<BinWord> l);



// !SECTION! Efficiently enumerating elements in a hyperplane 

struct DifferentialPair {
 BinWord x;
 BinWord y;
};

class DifferentialPairEnumerator
{
private:
    BinWord a;
    BinWord max_bit;
    BinWord i;

public:
    DifferentialPairEnumerator(BinWord size, BinWord a);
    DifferentialPair next();
    bool ended();
};

inline DifferentialPairEnumerator::DifferentialPairEnumerator(BinWord size, BinWord _a)
{
    a = _a;
    max_bit = size/2;
    i = 0;
}

inline bool DifferentialPairEnumerator::ended()
{
    return i >= max_bit;
}

inline DifferentialPair DifferentialPairEnumerator::next()
{
    DifferentialPair n;
    n.x = i;
    n.y = i^a;
    i++;
    if (n.y < n.x)
    {
        n.x+=max_bit;
        n.y+=max_bit;
    }
    return n;
}

#endif
