#ifndef _f2_funcs_
#define _f2_funcs_

#include "sboxU.hpp"


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

BinWord cpp_linear_combination (const std::vector<BinWord> & l, BinWord mask);

Integer cpp_rank_of_vector_set(std::vector<BinWord> l);


class cpp_Linear_basis
{
private:
    std::map<Integer, BinWord> basis;
public:
    cpp_Linear_basis() : basis() {} ;

    cpp_Linear_basis(const std::vector<BinWord> & l) ;

    void add_to_span(BinWord x);

    bool is_in_span(BinWord x) const ;

    std::vector<BinWord> get_basis() const;

    inline Integer rank() const
    {
        return basis.size();
    }
    
    std::vector<BinWord> span() const;
};
     
#endif
