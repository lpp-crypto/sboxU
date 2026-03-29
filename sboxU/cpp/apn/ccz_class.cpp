#include "./ccz_class.hpp"


// Returns all EA automorphisms of the quadratic APN function s.
// See ccz_class.hpp for a full description of the (L, delta) pair format and
// how to convert to the canonical (A_affine, B) representation for comparison
// with cpp_equivalences_from_lat / ea_equivalences.
std::vector<std::pair<cpp_F2AffineMap, BinWord>> cpp_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads
    )
{
    cpp_S_box o = cpp_ortho_derivative(s);
    // Each ortho-derivative automorphism encodes a candidate (L_A_inv, L_B_T)
    // pair for the EA automorphisms of s.  The LAT-based search gives
    // lower-triangular matrices; block_A = L_A_inv, block_B = L_B.
    std::vector<cpp_F2AffineMap> automorphisms_asd =
        cpp_equivalences_from_lat(o, o, false, n_threads, "linear");
    BinWord pw_n = s.input_space_size();
    std::vector<std::pair<cpp_F2AffineMap, BinWord>> automorphisms;
    for(auto autom : automorphisms_asd)
    {
        auto abcd = cpp_ccz_block_decomposition(autom);
        cpp_F2AffineMap
            L_A_inv = cpp_F2AffineMap(abcd[0]),
            L_B = cpp_F2AffineMap(abcd[1]),
            L_A = L_A_inv.inverse(),
            L_B_T = L_B.transpose();
        cpp_S_box
            L_A_inv_sb = L_A_inv.get_cpp_S_box(),
            L_B_T_sb   = L_B_T.get_cpp_S_box(),
            L_A_sb = L_A.get_cpp_S_box(),
            L_B_sb = L_B.get_cpp_S_box();

        // Sanity check: L_B and L_A are automorphisms of the ortho-derivative.
        // (This is now also checked inside cpp_equivalences_from_lat, but kept
        // here as an explicit guard.)
        if (L_B_sb * o * L_A_sb != o)
            std::cout << "[ERROR] automorphisms of the ortho-derivative are actually not automorphisms!" << std::endl;

        cpp_FunctionGraph G_s(s);

        // For each candidate input shift delta, derive L_C (the linear
        // correction term) and build the upper-triangular EA mapping
        //   L = cpp_EA_mapping(L_A_inv, L_B_T, L_C)
        //     = [ L_B_T  |  L_C * L_A ]
        //       [   0    |     L_A    ]
        // then accept it only if it maps graph(s) to itself (XOR-equiv check).
        // The accepted (L, delta) pair is returned so callers can reconstruct
        // the canonical (A_affine, B) form:
        //   A_affine(x) = L_A_inv(x) XOR delta
        //   B(x)        = L_B_T(x)
        for(BinWord delta=0; delta<pw_n; delta++)
        {
            // Build L_C: the unique linear map satisfying the upper-triangular
            // EA condition for this (L_A_inv, L_B_T, delta) triple.
            // C_0 is the affine constant; img_L_C[i] = L_C(e_i).
            std::vector<BinWord> img_L_C(s.get_input_length(), 0);
            BinWord C_0 = s[delta] ^ L_B_T_sb[s[L_A_inv_sb[0]]];
            for (unsigned int i=0; i<s.get_input_length(); i++)
            {
                BinWord e_i = 1 << i;
                img_L_C[i] = C_0 ^ s[e_i^delta] ^ L_B_T_sb[s[L_A_inv_sb[e_i]]];
            }
            cpp_F2AffineMap
                L_C(img_L_C),
                L = cpp_EA_mapping(L_A_inv, L_B_T, L_C);
            // Accept only if graph(s) maps to itself under L (up to XOR shift).
            std::vector<BinWord> offsets = G_s.xor_equivalence(G_s.image_by(L));
            if (offsets.size() > 0)
                automorphisms.push_back({L, delta});
        }
    }
    return automorphisms;
}


std::vector<cpp_F2AffineMap> cpp_ea_mappings_from_ortho_derivative(
    const cpp_S_box & s,
    const cpp_S_box & s_prime,
    const unsigned int n_threads
    )
{
    cpp_S_box
        o = cpp_ortho_derivative(s),
        o_prime = cpp_ortho_derivative(s_prime);
    std::vector<cpp_F2AffineMap> automorphisms_asd =
        cpp_equivalences_from_lat(o, o_prime, false, n_threads, "linear");
    BinWord pw_n = s.input_space_size();
    std::vector<cpp_F2AffineMap> automorphisms;
    for(auto autom : automorphisms_asd)
    {   
        auto abcd = cpp_ccz_block_decomposition(autom);
        cpp_F2AffineMap
            L_A_inv = cpp_F2AffineMap(abcd[0]),
            L_B = cpp_F2AffineMap(abcd[1]),
            L_A = L_A_inv.inverse(),
            L_B_T = L_B.transpose();
        cpp_S_box
            L_A_inv_sb = L_A_inv.get_cpp_S_box(),
            L_B_T_sb   = L_B_T.get_cpp_S_box(),
            L_A_sb = L_A.get_cpp_S_box(),
            L_B_sb = L_B.get_cpp_S_box();
        
        
        cpp_FunctionGraph
            G_s(s),
            G_s_prime(s_prime);

        // now need to find C
        for(BinWord delta=0; delta<pw_n; delta++)
        {
            // we compute the image vectors of the canonical basis
            std::vector<BinWord> img_L_C(s.get_input_length(), 0);
            BinWord C_0 = s[delta] ^ L_B_T_sb[s_prime[L_A_inv_sb[0]]];
            for (unsigned int i=0; i<s.get_input_length(); i++)
            {
                BinWord e_i = 1 << i;
                img_L_C[i] = C_0 ^ s[e_i^delta] ^ L_B_T_sb[s_prime[L_A_inv_sb[e_i]]];
            }
            // then we deduce what would be the linear part of an automorphism
            cpp_F2AffineMap
                L_C(img_L_C),
                L = cpp_EA_mapping(L_A_inv, L_B_T, L_C);
            // and then we check if the resulting graphs are XOR-equivalent
            std::vector<BinWord> offsets = G_s.xor_equivalence(G_s_prime.image_by(L));
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
    // Extract just the matrices from the (matrix, delta) pairs; init_mappings
    // only needs the matrices (the delta shifts are not used here).
    auto automs_with_delta = cpp_automorphisms_from_ortho_derivative(s, n_threads);
    std::vector<cpp_F2AffineMap> automs;
    automs.reserve(automs_with_delta.size());
    for (auto & [L, delta] : automs_with_delta)
        automs.push_back(L);
    ws.init_mappings(automs);
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
