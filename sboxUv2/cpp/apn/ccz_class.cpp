#include "./ccz_class.hpp"



std::vector<cpp_BinLinearMap> cpp_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads
    )
{
    cpp_S_box o = cpp_ortho_derivative(s);
    auto lat_ortho = cpp_lat(o);
    vector<vector<cpp_BinLinearMap>> autom_blocks =
        cpp_equivalences_from_lat(
            lat_ortho, lat_ortho, false, n_threads, "linear");
    BinWord pw_n = s.input_space_size();
    std::vector<cpp_BinLinearMap> automorphisms;
    for(auto abcd : autom_blocks)
    {
        cpp_BinLinearMap
            L_A_inv = cpp_BinLinearMap(abcd[0]).transpose(),
            L_B_T = cpp_BinLinearMap(abcd[1]),
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
        cpp_FunctionGraph G_s(s);

        // now need to find C
        for(BinWord delta=0; delta<pw_n; delta++)
        {
            cpp_S_box
                add_delta = cpp_translation(delta, s.get_input_length()),
                C = s*add_delta + L_B_T_sb * s * L_A_inv_sb;
            cpp_Spectrum diff = cpp_differential_spectrum(C, n_threads);
            if (diff.contains(pw_n) && (diff[pw_n] == (pw_n-1)))
            {
                cpp_S_box
                    C_0 = cpp_translation(C[0], s.get_input_length());                
                cpp_BinLinearMap
                    L_C(C + C_0),
                    L = cpp_EA_mapping(L_A_inv, L_B_T, L_C);
                // sanity check
                std::cout << std::hex << delta << "  "
                          << " "
                          << G_s.get_ccz_equivalent_function(L).content_string_repr()
                          << std::endl ;
                    
                // !CONTINUE! The automorphism for the ortho-derivative is correct, but not for the function
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
    
