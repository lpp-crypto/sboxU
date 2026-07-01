#include "./ea_from_vq.hpp"
#include <map>

std::map<cpp_BinLinearBasis, std::vector<unsigned int>> cpp_image_of_space_by_group(
    std::vector<cpp_F2AffineMap> G, cpp_BinLinearBasis V)
{
    std::map<cpp_BinLinearBasis, std::vector<unsigned int>> images;
    for (unsigned int i = 0; i < G.size(); i++)
        images[V.image_by(G[i])].push_back(i);
    return images;
}

/// @brief Check if any g1(V1) is equal to g2(V2) for any g1 in G1 and g2 in G2,returns corresponding elements if true
/// @param G1 a group of affine maps as a std::vector<cpp_F2AffineMap>
/// @param G2 a group of affine maps as a std::vector<cpp_F2AffineMap>
/// @param V1 a vector space as a cpp_BinLinearBasis
/// @param V2 a vector space as a cpp_BinLinearBasis
/// @return vector of cpp_F2AffineMap that satisfies L(V1) = V2
std::vector<cpp_F2AffineMap> cpp_product_walsh_match(
    std::vector<cpp_F2AffineMap> G1,
    std::vector<cpp_F2AffineMap> G2,
    cpp_BinLinearBasis V1,
    cpp_BinLinearBasis V2)
{
    std::vector<cpp_F2AffineMap> result;

    auto images1 = cpp_image_of_space_by_group(G1, V1);
    auto images2 = cpp_image_of_space_by_group(G2, V2);

    for (auto& [W, idx1] : images1) {
        auto it = images2.find(W);
        if (it == images2.end()) continue;
        for (unsigned int i : idx1)
            for (unsigned int j : it->second)
                result.push_back(G2[j].inverse() * G1[i]);
    }

    return result;
}

/// @brief Early-abort test: check if any pair (g1, g2) satisfies g1(V1) = g2(V2).
/// Builds the image map for the larger group (fewer stream iterations in expectation),
/// then streams the smaller group and returns true as soon as a matching image is found.
bool cpp_product_walsh_match_any(
    const std::vector<cpp_F2AffineMap>& G1,
    const std::vector<cpp_F2AffineMap>& G2,
    const cpp_BinLinearBasis& V1,
    const cpp_BinLinearBasis& V2)
{
    if (G1.size() >= G2.size()) {
        auto images = cpp_image_of_space_by_group(G1, V1);
        for (const auto& g2 : G2)
            if (images.count(V2.image_by(g2))) return true;
    } else {
        auto images = cpp_image_of_space_by_group(G2, V2);
        for (const auto& g1 : G1)
            if (images.count(V1.image_by(g1))) return true;
    }
    return false;
}


/// @brief Test whether f and g are EA-equivalent using Walsh zero spaces.
///
/// Implements the property: f ≃_EA g iff V_f and V_g lie in the same Aut(q_f)^T
/// orbit in WS(q_f), where q_f is the quadratic CCZ representative of f, and
/// V_f, V_g are the Walsh zero spaces of q_f corresponding to f and g respectively.
///
/// @param f    An S-box expected to be in the CCZ class of a quadratic APN function.
/// @param g    An S-box expected to be in the CCZ class of a quadratic APN function.
/// @param n_threads  Number of threads for parallel computation.
/// @param mode "standard": iterate over the full Aut(q_f).
///             "product":  use the G1 ⋊ G2 semidirect-product structure of Aut(q_f)
///                         via cpp_product_walsh_match.
/// @return A vector containing an EA mapping from q_f to q_g if f ≃_EA g, empty otherwise.
std::vector<cpp_F2AffineMap> cpp_ea_mapping_from_vq(
    const cpp_S_box f,
    const cpp_S_box g,
    const unsigned int n_threads,
    const std::string & mode)
{
    // --- Find the quadratic representative of f ---
    cpp_WalshZeroesSpaces WS_f(f, n_threads);
    WS_f.init_mappings();
    cpp_FunctionGraph graph_f(f);

    cpp_F2AffineMap map_q_f;
    cpp_S_box q_f;
    int idx_f = -1;
    for (int i = 0; i < (int)WS_f.mappings.size(); i++) {
        cpp_S_box temp = graph_f.get_ccz_equivalent_function(WS_f.mappings[i]);
        if (!cpp_is_degree_bigger_than(temp, 2)) {
            map_q_f = WS_f.mappings[i];
            q_f = temp;
            idx_f = i;
            break;
        }
    }
    if (idx_f < 0) return {};

    // --- Find the quadratic representative of g ---
    cpp_WalshZeroesSpaces WS_g(g, n_threads);
    WS_g.init_mappings();
    cpp_FunctionGraph graph_g(g);

    cpp_F2AffineMap map_q_g;
    cpp_S_box q_g;
    int idx_g = -1;
    for (int i = 0; i < (int)WS_g.mappings.size(); i++) {
        cpp_S_box temp = graph_g.get_ccz_equivalent_function(WS_g.mappings[i]);
        if (!cpp_is_degree_bigger_than(temp, 2)) {
            map_q_g = WS_g.mappings[i];
            q_g = temp;
            idx_g = i;
            break;
        }
    }
    if (idx_g < 0) return {};

    // --- Get an EA map q_f → q_g (empty if not CCZ-equivalent) ---
    auto ea_q = cpp_ea_mappings_from_ortho_derivative(q_f, q_g, n_threads);
    if (ea_q.empty()) return {};

    // --- Express both functions as Walsh zero spaces in WS(q_f) ---
    //
    // map_q_f maps Graph(f) → Graph(q_f), so L^{-T}_{map_q_f} maps WS_f → WS_{q_f}
    cpp_BinLinearBasis V_f = WS_f.bases[idx_f].image_by(map_q_f.inverse().transpose());

    // map_q_g maps Graph(g) → Graph(q_g), so L^{-T}_{map_q_g} maps WS_g → WS_{q_g}
    // ea_q[0] maps Graph(q_f) → Graph(q_g), so L^T_{ea_q[0]} maps WS_{q_g} → WS_{q_f}
    cpp_BinLinearBasis V_g = WS_g.bases[idx_g]
        .image_by(map_q_g.inverse().transpose())
        .image_by(ea_q[0].transpose());

    // --- f ≃_EA g iff V_f and V_g are in the same Aut(q_f)^T orbit ---
    if (mode == "product") {
        // Exploit the semidirect product Aut(q_f) = G1 ⋊ G2:
        //   G1 = EL graph automorphisms, G2 = derivative automorphisms
        auto G1 = cpp_graph_el_automorphisms_from_ortho_derivative(q_f, n_threads);
        auto G2 = cpp_graph_automorphisms_from_derivatives(q_f);
        // G1_t  = lin(g1)^T    for g1 in G1
        // G2_ti = lin(g2)^{-T} for g2 in G2
        std::vector<cpp_F2AffineMap> G1_t, G2_ti;
        for (auto & g1 : G1)
            G1_t.push_back((g1 + g1.get_cstte()).transpose());
        for (auto & g2 : G2)
            G2_ti.push_back((g2 + g2.get_cstte()).transpose().inverse());
        if (cpp_product_walsh_match_any(G1_t, G2_ti, V_f, V_g))
            return {ea_q[0]};
    } else {
        auto Aut_q = cpp_automorphisms_from_ortho_derivative(q_f, n_threads);
        for (auto & B : Aut_q)
            if (V_f.image_by(B.transpose()) == V_g)
                return {ea_q[0]};
    }
    return {};
}
