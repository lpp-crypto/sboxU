#include "./ccz_class.hpp"
#include "ortho_derivative.hpp"
#include "sn.hpp"
#include <set>


// !TODO! Unifomize the output :
//   - the current output of cpp_equivalences_from_lat is lower-triangular
//   - the output of cpp_automorphisms_from_ortho_derivative is upper-triangular
// !TODO! change this function into cpp_ea_equivalence_from_ortho_derivate ?
std::vector<cpp_F2AffineMap> cpp_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads,
    const std::string & mode
    )
{
    if (mode == "product") {
        std::vector<cpp_F2AffineMap> G1 = cpp_graph_el_automorphisms_from_ortho_derivative(s, n_threads);
        std::vector<cpp_F2AffineMap> G2 = cpp_graph_automorphisms_from_derivatives(s);
        std::vector<cpp_F2AffineMap> result;
        result.reserve(G1.size() * G2.size());
        for (const auto & g1 : G1)
            for (const auto & g2 : G2)
                result.push_back(g1 * g2);
        return result;
    }
    // Defaults to previous implementation if no mode is provided
    cpp_S_box o = cpp_ortho_derivative(s);
    std::vector<cpp_F2AffineMap> automorphisms_asd =
        cpp_equivalences_from_lat(o, o, false, n_threads, "linear");
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
            cpp_F2AffineMap
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



/// @brief Compute the trivial graph automorphisms that correspond to the derivatives of a quadratic functions.
/// @param s a cpp_S_box
/// @return a vector of ccp_F2AffineMpa that represents the automorphism group from the derivatives of s
std::vector<cpp_F2AffineMap>cpp_graph_automorphisms_from_derivatives(cpp_S_box s){

    std::vector<cpp_F2AffineMap> derivative_group;
    // We work with s_0, the function translated at 0 in 0
    cpp_S_box s_0;
    if(s[0]==0){
        s_0 = s;
    }
    else{
        s_0 =  cpp_translation(s[0],s.get_input_length())*s;
    }
    std::vector<BinWord> Zero_lut(s.get_lut().size(),0);
    cpp_F2AffineMap Id_n = identity_F2AffineMap(s.get_input_length());
    cpp_F2AffineMap Id_2n = identity_F2AffineMap(2*s.get_input_length());
    cpp_S_box derivative;
    cpp_F2AffineMap aut;
    // We use the shape A: x -> x+a, B: x -> x+F[a], C: x -> DaF[x] + s[a] and D = 0
    for (BinWord a = 0; a < s.get_lut().size(); a++) {
        cpp_S_box derivative = s_0.derivative(a);
        if (cpp_algebraic_degree(derivative) <= 1) {
            aut = cpp_F2AffineMap_from_blocks(
                Id_n+a,
                Id_n+s_0[a],
                cpp_F2AffineMap_from_lut(derivative)+s_0[a],
                cpp_F2AffineMap_from_lut(Zero_lut)
            );
            // We put our elements back to the correct function by conjugacy
            derivative_group.push_back((Id_2n+s[0])*aut*(Id_2n+s[0]));
        }
    }
    return derivative_group;
}


/// @brief Computes the EL graph automorphisms of a quadratic function using its ortho-derivative 
/// @param s a cpp_S_box
/// @param n_threads number of threads for parallel computation
/// @return a vector of cpp_F2AffineMap that represents the EL automorphism group of s
std::vector<cpp_F2AffineMap> cpp_graph_el_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads
    )
{
    std::vector<BinWord> Zero_lut(s.get_lut().size(),0);
    cpp_F2AffineMap 
        Zero = cpp_F2AffineMap_from_lut(Zero_lut),
        Id_2n = identity_F2AffineMap(2*s.get_input_length());

    cpp_S_box s_0;
    // We work with s_0, the function translated at 0 in 0
    if(s[0]==0){
        s_0 = s;
    }
    else{
        s_0 =  cpp_translation(s[0],s.get_input_length())*s;
    }

    cpp_S_box o = cpp_ortho_derivative(s_0);
    std::vector<cpp_F2AffineMap> 
        automorphisms_pi = cpp_equivalences_from_lat(o,o,false,n_threads,"linear"),
        automorphisms;

    cpp_F2AffineMap aut;
    for(auto autom : automorphisms_pi)
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

        // Careful with all the inverse-transpose, this formula should be good
        cpp_S_box C_candidate = s_0*L_A_sb + L_B_T_sb*s_0;
        if (cpp_algebraic_degree(C_candidate) <=1){
            aut = cpp_F2AffineMap_from_blocks(
                L_A,
                L_B_T,
                cpp_F2AffineMap_from_lut(C_candidate),
                Zero);
            // We put our elements back to the correct function by conjugacy
            automorphisms.push_back((Id_2n+s[0])*aut*(Id_2n+s[0]));
        }
    }
    return automorphisms;
}


/// @brief Computes a set of generators of the group of derivative automorphisms of an APN quadratic
/// function. Identical to cpp_graph_automorphisms_from_derivatives, except that the automorphism
/// associated to a direction a is only computed for a in the canonical basis (a = 1, 2, 4, ...).
/// Since a -> automorphism(a) is a group homomorphism from (F_2^n, XOR) to the derivative
/// automorphism group, and the canonical basis generates F_2^n, the n resulting automorphisms
/// generate the full group returned by cpp_graph_automorphisms_from_derivatives.
/// @param s a cpp_S_box
/// @return a vector of cpp_F2AffineMap, a generating set of the group of derivative automorphisms of s
std::vector<cpp_F2AffineMap> cpp_gen_set_graph_automorphisms_from_derivative(
    const cpp_S_box & s
    )
{
    std::vector<cpp_F2AffineMap> gen_set;
    // We work with s_0, the function translated at 0 in 0
    cpp_S_box s_0;
    if(s[0]==0){
        s_0 = s;
    }
    else{
        s_0 =  cpp_translation(s[0],s.get_input_length())*s;
    }
    std::vector<BinWord> Zero_lut(s.get_lut().size(),0);
    cpp_F2AffineMap Id_n = identity_F2AffineMap(s.get_input_length());
    cpp_F2AffineMap Id_2n = identity_F2AffineMap(2*s.get_input_length());
    cpp_F2AffineMap aut;
    // We use the shape A: x -> x+a, B: x -> x+F[a], C: x -> DaF[x] + s[a] and D = 0,
    // restricting a to the canonical basis e_0, ..., e_{n-1}
    for (unsigned int i = 0; i < s.get_input_length(); i++) {
        BinWord a = (BinWord)1 << i;
        cpp_S_box derivative = s_0.derivative(a);
        if (cpp_algebraic_degree(derivative) <= 1) {
            aut = cpp_F2AffineMap_from_blocks(
                Id_n+a,
                Id_n+s_0[a],
                cpp_F2AffineMap_from_lut(derivative)+s_0[a],
                cpp_F2AffineMap_from_lut(Zero_lut)
            );
            // We put our elements back to the correct function by conjugacy
            gen_set.push_back((Id_2n+s[0])*aut*(Id_2n+s[0]));
        }
    }
    return gen_set;
}


// !SECTION! Helpers for cpp_gen_set_F2AffineMap_group


/// @brief Canonical, orderable representation of an F2AffineMap (its image vectors plus its
/// constant term), used as a std::set key in the subgroup closure below.
static std::vector<BinWord> cpp_affine_key(const cpp_F2AffineMap & L)
{
    std::vector<BinWord> key = L.get_image_vectors();
    key.push_back(L.get_cstte());
    return key;
}


/// @brief Naive primality test by trial division (group orders here are small enough for this to be instant).
static bool cpp_is_prime(size_t m)
{
    if (m < 2) return false;
    if (m % 2 == 0) return m == 2;
    for (size_t d = 3; d*d <= m; d += 2)
        if (m % d == 0) return false;
    return true;
}


/// @brief Enumerates <gens> by closing {id} under right-multiplication by gens. In a finite group
/// this yields the full subgroup generated by gens (inverses come for free, since g^{-1} =
/// g^{ord(g)-1}). This is the orbit/BFS form of Dimino's algorithm (Butler 1991).
static std::set<std::vector<BinWord>> cpp_subgroup_closure(
    const std::vector<cpp_F2AffineMap> & gens,
    const cpp_F2AffineMap & id)
{
    std::set<std::vector<BinWord>> seen;
    seen.insert(cpp_affine_key(id));
    std::vector<cpp_F2AffineMap> frontier = {id};
    while (!frontier.empty())
    {
        std::vector<cpp_F2AffineMap> next_frontier;
        for (const auto & x : frontier)
            for (const auto & g : gens)
            {
                cpp_F2AffineMap y = x * g;
                if (seen.insert(cpp_affine_key(y)).second)
                    next_frontier.push_back(y);
            }
        frontier = next_frontier;
    }
    return seen;
}


/// @brief Checks whether <candidate> generates exactly G (as sets, via subgroup closure).
static bool cpp_is_generating_set(
    const std::vector<cpp_F2AffineMap> & candidate,
    const std::vector<cpp_F2AffineMap> & G,
    const cpp_F2AffineMap & id)
{
    std::set<std::vector<BinWord>> closure = cpp_subgroup_closure(candidate, id);
    if (closure.size() != G.size())
        return false;
    for (const auto & g : G)
        if (closure.find(cpp_affine_key(g)) == closure.end())
            return false;
    return true;
}


/// @brief Computes a set of generators of a group of F2AffineMap elements depending on its properties.
///
/// - |G| == 1: G only contains the identity, so it is returned as-is.
/// - |G| prime: any non-identity element generates G on its own.
/// - otherwise: draws floor(log2(|G|))+1 elements of G uniformly at random. In "probabilistic" mode
///   they are returned immediately (a generating set with high probability, unverified). In
///   "deterministic" mode (default) they are accepted only if their subgroup closure is exactly G;
///   otherwise a fresh batch is drawn and checked again.
///
/// @param G group represented by the list of its elements
/// @param mode "deterministic" (verify before returning) or "probabilistic" (return unverified)
/// @return a vector of cpp_F2AffineMap that represents a set of generators of the group
std::vector<cpp_F2AffineMap> cpp_gen_set_F2AffineMap_group(
    std::vector<cpp_F2AffineMap> G,
    const std::string & mode)
{
    std::vector<cpp_F2AffineMap> gen_set;
    if (G.empty())
        return gen_set;

    cpp_F2AffineMap id = identity_F2AffineMap(G[0].get_input_length());

    if (G.size() == 1)
        return G;

    if (cpp_is_prime(G.size()))
    {
        for (const auto & g : G)
            if (!(g == id))
            {
                gen_set.push_back(g);
                return gen_set;
            }
        return gen_set;
    }

    unsigned int n_draw = (unsigned int)std::log2((double)G.size()) + 1;
    while (true)
    {
        std::vector<cpp_F2AffineMap> candidate;
        for (unsigned int i = 0; i < n_draw; i++)
            candidate.push_back(G[switching_neighbors_rand_int_64_cpp() % G.size()]);

        if (mode == "probabilistic")
            return candidate;

        if (cpp_is_generating_set(candidate, G, id))
            return candidate;
    }
}




std::vector<cpp_S_box> cpp_enumerate_ea_classes_quadratic_apn(
    const cpp_S_box &s,
    const unsigned int n_threads,
    const std::string & mode
    )
{
    // initializing Walsh zeroes
    cpp_WalshZeroesSpaces ws(s, n_threads);
    // computing linear automorphisms and reducing Walsh zeroes spaces
    if (mode == "product")
        ws.init_mappings(
            cpp_graph_el_automorphisms_from_ortho_derivative(s, n_threads),
            cpp_graph_automorphisms_from_derivatives(s)
            );
    else if (mode == "generators")
    {
        // G1 = EL automorphisms, G2 = derivative automorphisms; their generating sets
        // concatenated generate Aut(q_f) = G1 ⋊ G2
        std::vector<cpp_F2AffineMap> gens = cpp_gen_set_F2AffineMap_group(
            cpp_graph_el_automorphisms_from_ortho_derivative(s, n_threads),
            "deterministic"
            );
        std::vector<cpp_F2AffineMap> gens2 = cpp_gen_set_graph_automorphisms_from_derivative(s);
        gens.insert(gens.end(), gens2.begin(), gens2.end());
        ws.init_mappings_generators(gens);
    }
    else
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
