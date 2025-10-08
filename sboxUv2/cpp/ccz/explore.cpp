#include "explore.hpp"


cpp_S_box cpp_ccz_equivalent_function(
    const cpp_S_box & s,
    const cpp_BinLinearMap & mapping
    )
{
    return cpp_FunctionGraph(s).get_ccz_equivalent_function(mapping);
}


std::vector<cpp_S_box> cpp_enumerate_ea_classes(
    const cpp_S_box & s,
    const unsigned int n_threads
    ) 
{
    cpp_WalshZeroesSpaces ws(s, n_threads);
    ws.init_mappings();
    return cpp_enumerate_ea_classes(s, ws);
}


std::vector<cpp_S_box> cpp_enumerate_ea_classes(
    const cpp_S_box & s,
    const cpp_WalshZeroesSpaces &ws
    ) 
{
    std::vector<cpp_S_box> result;
    cpp_FunctionGraph g(s);
    for(auto &L : ws.mappings)
        result.push_back(g.get_ccz_equivalent_function(L));
    return result;
}
