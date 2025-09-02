#include "zeroes.hpp"


cpp_WalshZeroesSpaces::cpp_WalshZeroesSpaces(
    const cpp_S_box & s,
    const unsigned int n_threads)
{
    // setting the machinery to parse BinWords into 2 vectors
    n = s.get_output_length();
    total_size = s.get_input_length() + s.get_output_length();
    mask = 0;
    for(unsigned int i=0; i<s.get_output_length(); i++)
        mask = (mask << 1) | 1;
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
    bases = cpp_extract_bases(
        zeroes_coord,
        n,
        n_threads,
        "fixed dimension"
        );
}


void cpp_WalshZeroesSpaces::init_mappings()
{
    // building the mappings by transposing
    for(auto &b : bases)
    {
        std::vector<BinWord> img = cpp_complete_basis(b, total_size);
        std::reverse(img.begin(), img.end());
        cpp_BinLinearMap L(img);
        mappings.push_back(L.transpose());
    }
}


void cpp_WalshZeroesSpaces::init_mappings(
    const std::vector<cpp_BinLinearMap> & automorphisms
    )
{
    // computing the image of each basis
    std::map<cpp_BinLinearBasis, unsigned int> preimages;
    for (unsigned int i=0; i<bases.size(); i++)
        preimages[cpp_BinLinearBasis(bases[i])] = i;
    // initializing walsh zeroes automorphisms
    std::vector<cpp_BinLinearMap> A;
    A.reserve(automorphisms.size());
    for(auto & a_i : automorphisms)
        A.push_back(a_i.transpose());
    // checking if an automorphism maps a space to another
    std::vector<bool> relevant(preimages.size(), true);
    for (auto & space : preimages)
        if (relevant[space.second])
            for (auto & Aj : A)
            {
                cpp_BinLinearBasis img = space.first.image_by(Aj);
                if (preimages.contains(img))
                {
                    unsigned int index = preimages[img];
                    if (index != space.second)
                        relevant[index] = false;
                }
            }
    // building the mappings by transposing
    for(unsigned int i=0; i<bases.size(); i++)
        if (relevant[i])
        {
            std::vector<BinWord> img = cpp_complete_basis(bases[i],
                                                          total_size);
            std::reverse(img.begin(), img.end());
            cpp_BinLinearMap L(img);
            mappings.push_back(L.transpose());
        }
}



cpp_Spectrum cpp_WalshZeroesSpaces::thickness_spectrum() const
{
    cpp_Spectrum result;
    for(auto &b : bases)
    {
        std::vector<BinWord> proj(b.size(), 0);
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
