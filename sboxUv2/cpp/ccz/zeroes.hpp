#ifndef _CCZ_ZEROES_
#define _CCZ_ZEROES_

#include "../sboxU.hpp"
#include "../core/include.hpp"
#include "../statistics/include.hpp"
#include "../algorithms/spaceSearch.hpp"



class cpp_WalshZeroesSpaces
{
private:
    BinWord mask;
    BinWord n;
    unsigned int total_size;
public:
    std::vector<std::vector<BinWord> > bases;
    std::vector<cpp_BinLinearMap> mappings;

    cpp_WalshZeroesSpaces(
        const cpp_S_box & s,
        const unsigned int n_threads
        );

    void init_mappings();

    void init_mappings(const std::vector<cpp_BinLinearMap> & automorphisms);
    
    cpp_Spectrum thickness_spectrum() const;
};


cpp_Spectrum cpp_thickness_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

#endif
