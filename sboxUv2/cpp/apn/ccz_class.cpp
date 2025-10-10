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
            L_A_inv = cpp_BinLinearMap(cpp_S_box(ab.first)).transpose(),
            L_B_T = cpp_BinLinearMap(cpp_S_box(ab.second)),
            L_A = L_A_inv.inverse(),
            L_B = L_B_T.transpose();
        cpp_S_box
            L_A_inv_sb = L_A_inv.get_cpp_S_box(),
            L_B_T_sb   = L_B_T.get_cpp_S_box(),
            L_A_sb = L_A.get_cpp_S_box(),
            L_B_sb = L_B.get_cpp_S_box();
        
        // sanity check
        if (L_B_sb*o*L_A_sb != o)
            std::cout << "[ERROR] automorphisms of the ortho-derivative are actually not automorphisms!" << std::endl;
        cpp_FunctionGraph G_o(o);

        // now need to find C
        for(BinWord delta=0; delta<s.input_space_size(); delta++)
        {
            cpp_S_box
                add_delta = cpp_translation(delta, s.get_input_length()),
                C = s*add_delta + L_B_T_sb * s * L_A_inv_sb;
            cpp_Spectrum diff = cpp_differential_spectrum(C, n_threads);
            if (diff[s.input_space_size()] == (s.input_space_size()-1))
            {
                cpp_S_box
                    C_0 = cpp_translation(C[0], s.get_input_length());                
                cpp_BinLinearMap
                    L_C(C + C_0),
                    L = cpp_EA_mapping(L_A_inv, L_B_T, L_C);
                // sanity check
                std::cout << std::hex << delta << "  "
                          << L.rank()
                          << " "
                          << (G_o.get_ccz_equivalent_function(L) == o)
                          << std::endl ;
                    
                // !CONTINUE! Add a test that the automorphism is indeed a graph automorphism
                automorphisms.push_back(L);
            }
        }
    }
    std::cout << "TOTAL: " << std::dec << automorphisms.size() << std::endl;
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
    
