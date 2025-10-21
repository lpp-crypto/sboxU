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


// !SECTION! Initiliazing the mappings

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
    std::cout << "reducing ["
              << std::dec << bases.size()
              << "] bases using ["
              << std::dec << automorphisms.size()
              << "] automorphisms" << std::endl ;
    // computing the image of each basis
    std::map<cpp_BinLinearBasis, unsigned int> preimages;
    for (unsigned int i=0; i<bases.size(); i++)
        preimages[bases[i]] = i;
    // initializing walsh zeroes automorphisms
    std::vector<cpp_BinLinearMap> A;
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
            cpp_BinLinearMap L(img);
            mappings.push_back(L.transpose());
        }
}

// !SECTION! Applying a linear permutation

cpp_WalshZeroesSpaces cpp_WalshZeroesSpaces::image_by(
    const cpp_BinLinearMap & L
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
