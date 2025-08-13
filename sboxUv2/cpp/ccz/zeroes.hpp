#ifndef _CCZ_ZEROES_
#define _CCZ_ZEROES_

#include "../sboxU.hpp"
#include "../core/include.hpp"
#include "../statistics/include.hpp"
#include "../algorithms/spaceSearch.hpp"



class cpp_WalshZeroesSpaces
{
private:
    std::vector<std::vector<BinWord> > bases;
    BinWord mask;
    BinWord n;

public:
    cpp_WalshZeroesSpaces(
        const cpp_S_box & s,
        const unsigned int n_threads
        );
    cpp_Spectrum thickness_spectrum();
};

cpp_Spectrum cpp_thickness_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

#endif
