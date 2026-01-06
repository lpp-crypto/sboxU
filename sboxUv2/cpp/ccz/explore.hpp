#ifndef _CCZ_EXPLORE_
#define _CCZ_EXPLORE_

#include "../sboxU.hpp"
#include "../core/include.hpp"
#include "./zeroes.hpp"
#include "./graph.hpp"


cpp_S_box cpp_ccz_equivalent_function(
    const cpp_S_box & s,
    const cpp_BinLinearMap & L
    );


std::vector<cpp_S_box> cpp_enumerate_ea_classes(
    const cpp_S_box & s,
    const unsigned int n_threads
    );


std::vector<cpp_S_box> cpp_enumerate_ea_classes(
    const cpp_S_box & s,
    const cpp_WalshZeroesSpaces &ws
    );

#endif
