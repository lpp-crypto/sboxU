#include "explore.hpp"


// !SECTION! Basic functions


cpp_S_box cpp_ccz_equivalent_function(
    const cpp_S_box & s,
    const cpp_F2AffineMap & mapping
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


// !SECTION! CCZ-equivalence to a permutation


// !SUBSECTION! The general case

std::vector<cpp_S_box> cpp_enumerate_permutations_in_ccz_class(
    const cpp_S_box &s,
    const cpp_WalshZeroesSpaces &ws
    )
{
    std::vector<cpp_S_box> result;
    cpp_FunctionGraph g(s);
    for(auto &b1 : ws.bases)
        for (auto &b2 : ws.bases)
            if (cpp_is_sum_full_rank(b1, b2))
            {
                std::vector<BinWord>
                    first_img = b1.get_basis(),
                    second_img= b2.get_basis();
                first_img.insert(first_img.end(), second_img.begin(), second_img.end());
                cpp_F2AffineMap L = cpp_F2AffineMap(first_img).transpose();
                if (L.rank() == (s.get_input_length() + s.get_output_length()))
                    result.push_back(g.get_ccz_equivalent_function(L));
            }
    return result;
}


std::vector<cpp_S_box> cpp_enumerate_permutations_in_ccz_class(
    const cpp_S_box & s,
    const unsigned int n_threads
    )
{
    cpp_WalshZeroesSpaces ws(s, n_threads);
    ws.init_mappings();
    return cpp_enumerate_permutations_in_ccz_class(s, ws);
}



// !SUBSECTION! The particular case of EA-equivalence


// !TODO! implement functions for EA-equivalence to a permutation 

// !TODO! implement finding permutations in EA-equivalence class
