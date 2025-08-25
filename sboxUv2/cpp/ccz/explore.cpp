#include "explore.hpp"


cpp_S_box cpp_ccz_equivalent_function(
    const cpp_S_box & s,
    const std::vector<BinWord> mapping
    )
{
    return cpp_FunctionGraph(s).apply_basis_change(mapping);
}


std::vector<cpp_S_box> cpp_enumerate_ea_classes(
    const cpp_S_box & s,
    const unsigned int n_threads
    ) 
{
    std::vector<cpp_S_box> result;
    cpp_FunctionGraph g(s);
    cpp_WalshZeroesSpaces ws(s, n_threads);
    ws.init_mappings();
    for(auto &L : ws.mappings)
        result.push_back(g.apply_basis_change(L));
    return result;
}
