#include "./ccz_class.hpp"



std::vector<cpp_BinLinearMap> cpp_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads
    )
{
    cpp_S_box o = cpp_ortho_derivative(s);
    vector<pair<vector<BinWord>, vector<BinWord>>> autom_luts =
        cpp_linear_automorphisms_from_lat(
            cpp_lat(o),
            "alt_partition_diag_mappings",
            n_threads
            );
    std::vector<cpp_BinLinearMap> automorphisms;
    for(auto ab : autom_luts)
    {
        cpp_BinLinearMap
            L_A = cpp_BinLinearMap(cpp_S_box(ab.first)).transpose(),
            L_B = cpp_BinLinearMap(cpp_S_box(ab.second));
        // now need to find C
        // sanity check [PASSED]
        L_A = L_A.inverse();
        L_B = L_B.transpose();
        for(BinWord x=0; x<o.output_space_size(); x++)
            if (L_B(o[L_A(x)]) != o[x])
                std::cout << "well shit" << std::endl;

        // !CONTINUE!  
    }
    return automorphisms;
}


std::vector<cpp_S_box> cpp_enumerate_ea_classes_quadratic_apn(
    const cpp_S_box &s,
    const unsigned int n_threads
    )
{
    // initializing Walsh zeroes
    cpp_WalshZeroesSpaces ws(s, n_threads);
    // computing linear automorphisms of the ortho-derivative
    ws.init_mappings(cpp_automorphisms_from_ortho_derivative(s, n_threads));
    return cpp_enumerate_ea_classes(s, ws);
}
    
