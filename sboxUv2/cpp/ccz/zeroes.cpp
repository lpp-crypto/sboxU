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


void init_mappings(std::vector<cpp_BinLinearMap> automorphisms)
{

    // !TODO! implement init_mappings. It needs:
    // ! - a method to apply a cpp_BinLinearMap to a cpp_LinearBasis
    // ! - a comparison operator for cpp_LinearBasis so they can be used as std::map keys
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
