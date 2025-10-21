#ifndef _CCZ_ZEROES_
#define _CCZ_ZEROES_

#include "../sboxU.hpp"
#include "../core/include.hpp"
#include "../statistics/include.hpp"
#include "../algorithms/include.hpp"



class cpp_WalshZeroesSpaces
{
private:
    BinWord mask;
    BinWord n;
    unsigned int total_size;
public:
    std::vector<cpp_BinLinearBasis> bases;
    std::vector<cpp_BinLinearMap> mappings;

    cpp_WalshZeroesSpaces() : bases(0), mappings(0) {} ;
    
    cpp_WalshZeroesSpaces(
        const std::vector<cpp_BinLinearBasis> & _bases,
        const unsigned int _n,
        const unsigned int _total_size
        ):
        mask((1 << _n) - 1),
        n(_n),
        total_size(_total_size),
        bases(_bases.begin(), _bases.end()),
        mappings(0)
    {};
    
    cpp_WalshZeroesSpaces(
        const cpp_S_box & s,
        const unsigned int n_threads
        );

    void init_mappings();

    void init_mappings(const std::vector<cpp_BinLinearMap> & automorphisms);

    cpp_WalshZeroesSpaces image_by(const cpp_BinLinearMap & L) const;
    
    cpp_Spectrum thickness_spectrum() const;
};


cpp_Spectrum cpp_thickness_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

#endif
