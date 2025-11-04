#include "./ccz_class.hpp"


// !TODO! Unifomize the output :
//   - the current output of cpp_equivalences_from_lat is lower-triangular
//   - the output of cpp_automorphisms_from_ortho_derivative is upper-triangular
// !TODO! change this function into cpp_ea_equivalence_from_ortho_derivate ?
std::vector<cpp_BinLinearMap> cpp_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads
    )
{
    cpp_S_box o = cpp_ortho_derivative(s);
    std::vector<cpp_BinLinearMap> automorphisms_asd =
        cpp_equivalences_from_lat(o, o, false, n_threads, "linear");
    BinWord pw_n = s.input_space_size();
    std::vector<cpp_BinLinearMap> automorphisms;
    for(auto autom : automorphisms_asd)
    {   
        auto abcd = cpp_ccz_block_decomposition(autom);
        cpp_BinLinearMap
            L_A_inv = cpp_BinLinearMap(abcd[0]),
            L_B = cpp_BinLinearMap(abcd[1]),
            L_A = L_A_inv.inverse(),
            L_B_T = L_B.transpose();
        cpp_S_box
            L_A_inv_sb = L_A_inv.get_cpp_S_box(),
            L_B_T_sb   = L_B_T.get_cpp_S_box(),
            L_A_sb = L_A.get_cpp_S_box(),
            L_B_sb = L_B.get_cpp_S_box();
        
        // sanity check
        // !TODO! Remove : the necessary check is now done in cpp_equivalences_from_lat
        if (L_B_sb * o * L_A_sb != o)
            std::cout << "[ERROR] automorphisms of the ortho-derivative are actually not automorphisms!" << std::endl;
        
        cpp_FunctionGraph G_s(s);

        // now need to find C
        for(BinWord delta=0; delta<pw_n; delta++)
        {
            // we compute the image vectors of the canonical basis
            std::vector<BinWord> img_L_C(s.get_input_length(), 0);
            BinWord C_0 = s[delta] ^ L_B_T_sb[s[L_A_inv_sb[0]]];
            for (unsigned int i=0; i<s.get_input_length(); i++)
            {
                BinWord e_i = 1 << i;
                img_L_C[i] = C_0 ^ s[e_i^delta] ^ L_B_T_sb[s[L_A_inv_sb[e_i]]];
            }
            // then we deduce what would be the linear part of an automorphism
            cpp_BinLinearMap
                L_C(img_L_C),
                L = cpp_EA_mapping(L_A_inv, L_B_T, L_C);
            // and then we check if the resulting graphs are XOR-equivalent
            std::vector<BinWord> offsets = G_s.xor_equivalence(G_s.image_by(L));
            if (offsets.size() > 0)
                automorphisms.push_back(L);
        }
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
    

cpp_S_box cpp_ccz_equivalent_quadratic_function(
    const cpp_S_box & s,
    const unsigned int n_threads
    )
{
    for (auto &f : cpp_enumerate_ea_classes(s, n_threads))
        if (cpp_algebraic_degree(f) == 2)
            return f;
    return cpp_empty_S_box();
}
