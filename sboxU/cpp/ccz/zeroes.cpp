#include "zeroes.hpp"


// !SECTION! Construction

cpp_WalshZeroesSpaces::cpp_WalshZeroesSpaces(
    const cpp_S_box & s,
    const unsigned int n_threads) :
    mappings(0)
{
    // setting the machinery to parse BinWords into 2 vectors
    n = s.get_output_length();
    total_size = s.get_input_length() + s.get_output_length();
    mask = (1 << s.get_output_length()) - 1;
    // going over the zeroes in the LAT, column by column
    std::vector<BinWord> zeroes_coord;
    for (unsigned int b=0; b<s.output_space_size(); b++)
    {
        std::vector<Integer> w = cpp_walsh_transform(s.component(b));
        for (unsigned int a=0; a<w.size(); a++)
            if (w[a] == 0)
                zeroes_coord.push_back((a << n) | b);
    }
    // finding vector spaces of dimension n
    bases = cpp_extract_BinLinearBases(
        zeroes_coord,
        n,
        n_threads,
        "fixed dimension"
        );
}


// !SECTION! Initializing the mappings

void cpp_WalshZeroesSpaces::init_mappings()
{
    // building the mappings by transposing
    for(auto &b : bases)
    {
        std::vector<BinWord> img = cpp_complete_basis(b, total_size);
        std::reverse(img.begin(), img.end());
        cpp_F2AffineMap L(img);
        mappings.push_back(L.transpose());
    }
}


void cpp_WalshZeroesSpaces::init_mappings(
    const std::vector<cpp_F2AffineMap> & automorphisms
    )
{
    // computing the image of each basis
    std::map<cpp_BinLinearBasis, unsigned int> preimages;
    for (unsigned int i=0; i<bases.size(); i++)
        preimages[bases[i]] = i;
    // initializing walsh zeroes automorphisms
    std::vector<cpp_F2AffineMap> A;
    A.reserve(automorphisms.size());
    for(auto & a_i : automorphisms)
        A.push_back(a_i.transpose());
    // checking if an automorphism maps a space to another
    std::vector<bool> relevant(bases.size(), true);
    for (auto & space_and_i : preimages)
    {
        cpp_BinLinearBasis space = space_and_i.first;
        unsigned int i = space_and_i.second;
        if (relevant[i])
            for (auto & Aj : A)
            {
                cpp_BinLinearBasis img = space.image_by(Aj);
                if (preimages.contains(img))
                {
                    unsigned int index = preimages[img];
                    if (index != i)
                        relevant[index] = false;
                }
            }
    }
    // building the mappings by transposing
    for(unsigned int i=0; i<bases.size(); i++)
        if (relevant[i])
        {
            std::vector<BinWord> img = cpp_complete_basis(bases[i],
                                                          total_size);
            std::reverse(img.begin(), img.end());
            cpp_F2AffineMap L(img);
            mappings.push_back(L.transpose());
        }
}


void cpp_WalshZeroesSpaces::init_mappings(
    const std::vector<cpp_F2AffineMap> & automorphisms_1,
    const std::vector<cpp_F2AffineMap> & automorphisms_2
    )
{
    // G1_t  = lin(g1).T          for g1 in G1
    // G2_ti = lin(g2).T^{-1}     for g2 in G2
    std::vector<cpp_F2AffineMap> G1_t, G2_ti;
    for (auto & a : automorphisms_1)
        G1_t.push_back((a + a.get_cstte()).transpose());
    for (auto & b : automorphisms_2)
        G2_ti.push_back((b + b.get_cstte()).transpose().inverse());

    // Use the smaller group for step 1 to minimise the preimage map size.
    const auto & G_map  = (G1_t.size() <= G2_ti.size()) ? G1_t : G2_ti;
    const auto & G_walk = (G1_t.size() <= G2_ti.size()) ? G2_ti : G1_t;

    // Build the preimage map: for each V_i and each g in G_map,
    // insert W = g(V_i) → i.  The first space to claim a given W keeps it;
    // any later claimant is G_map-equivalent to the owner and is marked non-relevant.
    std::map<cpp_BinLinearBasis, unsigned int> preimage;
    std::vector<bool> relevant(bases.size(), true);
    for (unsigned int i = 0; i < bases.size(); i++) {
        for (auto & g : G_map) {
            auto [it, inserted] = preimage.emplace(bases[i].image_by(g), i);
            if (!inserted && it->second != i)
                relevant[i] = false;
        }
    }

    // For each relevant space, apply every g in G_walk to its preimage entries;
    // if the result is another preimage, mark that space non-relevant.
    for (auto & [W, i] : preimage) {
        if (!relevant[i]) continue;
        for (auto & g : G_walk) {
            auto it = preimage.find(W.image_by(g));
            if (it != preimage.end() && it->second != i)
                relevant[it->second] = false;
        }
    }

    // Build one mapping per relevant space
    for (unsigned int i = 0; i < bases.size(); i++)
        if (relevant[i])
        {
            std::vector<BinWord> img = cpp_complete_basis(bases[i], total_size);
            std::reverse(img.begin(), img.end());
            cpp_F2AffineMap L(img);
            mappings.push_back(L.transpose());
        }
}



// !SECTION! Orbits under a group action

/// @brief Computes the orbits of the Walsh zero spaces of `ws` under the group generated by
/// `generators`, the action of a generator g on a space V being V.image_by(g).
///
/// For every (generator, space) pair, the image of the space is looked up among the original
/// spaces; if it is one of them, the two spaces are merged into the same orbit. This is done
/// with a union-find structure (path halving + union by size), so the whole computation is a
/// single pass over generators x spaces, each step doing one image_by call and one O(log N)
/// map lookup, with the orbits themselves recovered in linear time at the end.
///
/// @param ws a cpp_WalshZeroesSpaces instance; ws.bases are the spaces to partition into orbits
/// @param generators a generating set of the acting group, as a vector of cpp_F2AffineMap
/// @return the orbits, each orbit being a vector of indices into ws.bases
std::vector<std::vector<unsigned int>> cpp_Walsh_zero_orbits(
    const cpp_WalshZeroesSpaces & ws,
    const std::vector<cpp_F2AffineMap> & generators)
{
    const std::vector<cpp_BinLinearBasis> & spaces = ws.bases;
    const unsigned int N = spaces.size();

    // index lookup: space -> its index in `spaces`
    std::map<cpp_BinLinearBasis, unsigned int> index_of;
    for (unsigned int i = 0; i < N; i++)
        index_of[spaces[i]] = i;

    // union-find over [0, N), path halving + union by size
    std::vector<unsigned int> parent(N), rank_size(N, 1);
    std::iota(parent.begin(), parent.end(), 0);

    auto find = [&](unsigned int x) -> unsigned int
    {
        while (parent[x] != x)
        {
            parent[x] = parent[parent[x]];   // path halving
            x = parent[x];
        }
        return x;
    };

    auto unite = [&](unsigned int a, unsigned int b)
    {
        a = find(a);
        b = find(b);
        if (a == b)
            return;
        if (rank_size[a] < rank_size[b])
            std::swap(a, b);
        parent[b] = a;
        rank_size[a] += rank_size[b];
    };

    // note the action of every generator on every space, unioning matches as we go
    for (auto & g : generators)
        for (unsigned int i = 0; i < N; i++)
        {
            auto it = index_of.find(spaces[i].image_by(g));
            if (it != index_of.end())
                unite(i, it->second);
        }

    // bucket the indices by orbit root
    std::vector<std::vector<unsigned int>> buckets(N);
    for (unsigned int i = 0; i < N; i++)
        buckets[find(i)].push_back(i);

    std::vector<std::vector<unsigned int>> orbits;
    for (auto & b : buckets)
        if (!b.empty())
            orbits.push_back(std::move(b));

    return orbits;
}


/// @brief Builds one mapping per orbit of cpp_Walsh_zero_orbits(*this, generators), each mapping
/// derived from an (arbitrary) representative of its orbit -- since the orbits already account
/// for the action of every generator, only one representative per orbit is needed instead of the
/// full per-space relevance scan that the other init_mappings overloads perform.
void cpp_WalshZeroesSpaces::init_mappings_generators(
    const std::vector<cpp_F2AffineMap> & generators
    )
{
    // Walsh zero spaces automorphisms act via the transpose; strip any constant term first
    // (transpose() throws on a non-zero constant), as generators of derivative/EL automorphism
    // groups may be affine rather than purely linear.
    std::vector<cpp_F2AffineMap> A;
    A.reserve(generators.size());
    for (auto & g : generators)
        A.push_back((g + g.get_cstte()).transpose());

    for (auto & orbit : cpp_Walsh_zero_orbits(*this, A))
    {
        std::vector<BinWord> img = cpp_complete_basis(bases[orbit[0]], total_size);
        std::reverse(img.begin(), img.end());
        cpp_F2AffineMap L(img);
        mappings.push_back(L.transpose());
    }
}



// !SECTION! Applying a linear permutation

cpp_WalshZeroesSpaces cpp_WalshZeroesSpaces::image_by(
    const cpp_F2AffineMap & L
    ) const
{
    std::vector<cpp_BinLinearBasis> img_bases;
    img_bases.reserve(bases.size());
    for(auto & b : bases)
        img_bases.push_back(b.image_by(L));
    return cpp_WalshZeroesSpaces(img_bases, n, total_size);
}



// !SECTION! Thickness spectrum 

cpp_Spectrum cpp_WalshZeroesSpaces::thickness_spectrum() const
{
    cpp_Spectrum result;
    for(auto &bl : bases)
    {
        std::vector<BinWord>
            b = bl.get_basis(),
            proj(b.size(), 0);
        for (unsigned int i=0; i<b.size(); i++)
            proj[i] = b[i] & mask;
        result.incr(cpp_rank_of_vector_set(proj));
    }
    return result;
}



cpp_Spectrum cpp_thickness_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads)
{
    cpp_WalshZeroesSpaces wz(s, n_threads);
    return wz.thickness_spectrum();
}
