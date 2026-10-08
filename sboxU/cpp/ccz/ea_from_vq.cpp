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
/// then iterates over the smaller group and returns true as soon as a matching image is found.
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


/// @brief Find (i, j) s.t. (G1[i]*G2[j])^{-T}(Vf) == Vg.
/// Equivalent condition: Vf.image_by(lin(G2[j])^{-T}) == Vg.image_by(lin(G1[i])^T).
/// Builds hashmap over the larger group; iterates over the smaller.
/// @return {i, j} or {-1, -1} if no matching pair exists.
std::pair<int,int> cpp_product_walsh_match_indices(
    const std::vector<cpp_F2AffineMap>& G1,
    const std::vector<cpp_F2AffineMap>& G2,
    const cpp_BinLinearBasis& Vf,
    const cpp_BinLinearBasis& Vg)
{
    auto lin_t = [](const cpp_F2AffineMap& g){ return (g + g.get_cstte()).transpose(); };
    if (G1.size() >= G2.size()) {
        std::map<cpp_BinLinearBasis, int> img1;
        for (int i = 0; i < (int)G1.size(); i++)
            img1[Vg.image_by(lin_t(G1[i]))] = i;
        for (int j = 0; j < (int)G2.size(); j++) {
            auto it = img1.find(Vf.image_by(lin_t(G2[j])));
            if (it != img1.end()) return {it->second, j};
        }
    } else {
        std::map<cpp_BinLinearBasis, int> img2;
        for (int j = 0; j < (int)G2.size(); j++)
            img2[Vf.image_by(lin_t(G2[j]))] = j;
        for (int i = 0; i < (int)G1.size(); i++) {
            auto it = img2.find(Vg.image_by(lin_t(G1[i])));
            if (it != img2.end()) return {i, it->second};
        }
    }
    return {-1, -1};
}

/// @brief Test whether f and g are EA-equivalent using Walsh zero spaces.
///
/// Implements the property: f =_EA g iff V_f and V_g lie in the same Aut(q_f)^T
/// orbit in WS(q_f), where q_f is the quadratic CCZ representative of f, and
/// V_f, V_g are the Walsh zero spaces of q_f corresponding to f and g respectively.
///
/// @param f    An S-box expected to be in the CCZ class of a quadratic APN function.
/// @param g    An S-box expected to be in the CCZ class of a quadratic APN function.
/// @param n_threads  Number of threads for parallel computation.
/// @param mode "standard": Iterate over the full Aut(q_f).
///             "product":  Use the G1 ⋊ G2 semidirect-product structure of Aut(q_f)
///                         via cpp_product_walsh_match.
///             "generators": Walk of the orbit of Vf under generating sets of
///                         G1 and G2, never enumerating Aut(q_f) itself.
/// @return A vector containing an EA mapping from q_f to q_g if f =_EA g, empty otherwise.
std::vector<cpp_F2AffineMap> cpp_ea_mapping_from_vq(
    const cpp_S_box f,
    const cpp_S_box g,
    const unsigned int n_threads,
    const std::string & mode)
{
    //Find the quadratic representative of f
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

    //Find the quadratic representative of g
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

    //f =_EA g iff Vf and Vg lie in the same Aut(q_f)^{-T} orbit
    if (mode == "product") {
        // Use the semidirect product Aut(q_f) = G1 ⋊ G2 to find (G1[i], G2[j])
        auto ea_maps = cpp_ea_mappings_from_ortho_derivative(q_g, q_f, n_threads, "single");
        if (ea_maps.empty()) return {};
        auto ea = ea_maps[0];
        auto G1 = cpp_graph_el_automorphisms_from_ortho_derivative(q_f, n_threads);
        auto G2 = cpp_graph_automorphisms_from_derivatives(q_f);
        auto mit_f = map_q_f.inverse().transpose();
        auto mit_g = map_q_g.inverse().transpose();
        cpp_BinLinearBasis Vf = WS_f.bases[idx_f].image_by(mit_f).image_by(mit_f);
        cpp_BinLinearBasis Vg = WS_g.bases[idx_g].image_by(mit_g).image_by(mit_g)
                                                  .image_by(ea.transpose());
        auto [i, j] = cpp_product_walsh_match_indices(G1, G2, Vf, Vg);
        if (i < 0) return {};
        return {map_q_g.inverse() * ea * G1[i] * G2[j] * map_q_f};
    } else if (mode == "generators") {
        auto ea_maps = cpp_ea_mappings_from_ortho_derivative(q_g, q_f, n_threads, "single");
        if (ea_maps.empty()) return {};
        auto ea = ea_maps[0];
        // Aut(q_f) = G1 ⋊ G2 is generated by a generating set of G1 together with the n
        // derivative automorphisms of the canonical basis directions.
        std::vector<cpp_F2AffineMap> gens = cpp_gen_set_F2AffineMap_group(
            cpp_graph_el_automorphisms_from_ortho_derivative(q_f, n_threads), "deterministic");
        std::vector<cpp_F2AffineMap> gens2 = cpp_gen_set_graph_automorphisms_from_derivative(q_f);
        gens.insert(gens.end(), gens2.begin(), gens2.end());
        std::vector<cpp_F2AffineMap> gens_t;
        for (auto & h : gens)
            gens_t.push_back((h + h.get_cstte()).transpose());

        auto mit_f = map_q_f.inverse().transpose();
        auto mit_g = map_q_g.inverse().transpose();
        cpp_BinLinearBasis Vf = WS_f.bases[idx_f].image_by(mit_f).image_by(mit_f);
        cpp_BinLinearBasis Vg = WS_g.bases[idx_g].image_by(mit_g).image_by(mit_g)
                                                  .image_by(ea.transpose());
        if (Vf == Vg)
            return {map_q_g.inverse() * ea * map_q_f};

        // Walk of the orbit of Vf, keeping for every reached space W an
        // automorphism M of q_f with lin(M)^T(Vf) == W. As with "standard", M^T(Vf) == Vg means
        // L = M^{-1} satisfies L^{-T}(Vf) == Vg, so only the final M is inverted.
        std::map<cpp_BinLinearBasis, cpp_F2AffineMap> reached;
        reached.emplace(Vf, identity_F2AffineMap(2 * f.get_input_length()));
        std::vector<cpp_BinLinearBasis> frontier = {Vf};
        while (!frontier.empty()) {
            std::vector<cpp_BinLinearBasis> next;
            for (auto & V : frontier) {
                const cpp_F2AffineMap M = reached.at(V);
                for (unsigned int k = 0; k < gens.size(); k++) {
                    cpp_BinLinearBasis W = V.image_by(gens_t[k]);
                    if (reached.contains(W))
                        continue;
                    cpp_F2AffineMap MW = M * gens[k];
                    if (W == Vg)
                        return {map_q_g.inverse() * ea * MW.inverse() * map_q_f};
                    reached.emplace(W, MW);
                    next.push_back(W);
                }
            }
            frontier = std::move(next);
        }
    } else {
        // Map every Walsh zero space to the one of q_f, then scan the full Aut(q_f).
        auto ea_maps = cpp_ea_mappings_from_ortho_derivative(q_g, q_f, n_threads, "single");
        if (ea_maps.empty()) return {};
        auto ea = ea_maps[0];
        auto Aut_q = cpp_automorphisms_from_ortho_derivative(q_f, n_threads);
        auto mit_f = map_q_f.inverse().transpose();
        auto mit_g = map_q_g.inverse().transpose();
        cpp_BinLinearBasis Vf = WS_f.bases[idx_f].image_by(mit_f).image_by(mit_f);
        cpp_BinLinearBasis Vg = WS_g.bases[idx_g].image_by(mit_g).image_by(mit_g)
                                                  .image_by(ea.transpose());
        // Aut_q is a group: L^{-T}(Vf) == Vg for some L iff M^T(Vf) == Vg for M = L^{-1}, so
        // only the matching automorphism gets inverted (inverse() is a brute-force search, does not work after n=11 for now).
        for (const auto& M : Aut_q)
            if (Vf.image_by(M.transpose()) == Vg)
                return {map_q_g.inverse() * ea * M.inverse() * map_q_f};
    }
    return {};
}
