#include "./f2affinemap.hpp"


cpp_F2AffineMap::cpp_F2AffineMap(const cpp_S_box & lut) :
    image_vectors(0),
    input_length(0),
    output_length(0) 
{
    cstte = lut[0];
    for(unsigned int pos=1; pos<lut.size(); pos<<=1)
    {
        BinWord
            y = lut[pos] ^ cstte,
            m = cpp_msb(y);
        image_vectors.push_back(y);
        input_length += 1;
        if (m > output_length)
            output_length = m;
    }
    output_length ++;
}


cpp_F2AffineMap cpp_F2AffineMap::operator*(const cpp_F2AffineMap & l) const
{
    if (l.get_output_length() != input_length)
        throw std::runtime_error("mismatched input/output length when multiplying F2AffineMap:s");
    else
    {
        std::vector<BinWord> images(l.get_input_length(), 0);
        for(unsigned int i=0; i<l.get_input_length(); i++)
            images[i] = cpp_linear_combination(image_vectors, l.image_vectors[i]);
        return cpp_F2AffineMap(
            images,
            l.get_input_length(),
            output_length,
            operator()(l.cstte)
            );
    }
};


cpp_F2AffineMap cpp_F2AffineMap::operator+(const cpp_F2AffineMap & l) const
{
    if (l.get_input_length() != input_length)
        throw std::runtime_error("mismatched input length when adding F2AffineMap:s");
    else
    {
        std::vector<BinWord> images(image_vectors.begin(),
                                    image_vectors.end());
        for(unsigned int i=0; i<images.size(); i++)
            images[i] ^= l.image_vectors[i];
        return cpp_F2AffineMap(
            images,
            input_length,
            output_length,
            cstte ^ l.cstte
            );
    }
}


cpp_F2AffineMap cpp_F2AffineMap::operator+(const BinWord cst) const
{
    return cpp_F2AffineMap(
        image_vectors,
        input_length,
        output_length,
        cst^cstte
        );
}


cpp_F2AffineMap cpp_F2AffineMap::operator|(const cpp_F2AffineMap & l) const
{
    
    std::vector<BinWord> images(image_vectors);
    for(auto v : l.image_vectors)
        images.push_back(v << output_length);
    return cpp_F2AffineMap(
        images,
        input_length + l.get_input_length(),
        cstte | (l.cstte << output_length)
        );
}


cpp_F2AffineMap cpp_F2AffineMap::inverse() const
{
    if (input_length <= 20)      // up to n=8, n^2.8 is greater than 2^n
                                 // but for now, we don't have anything else
    {
        // bruteforcing inputs and sorting those mapped to
        // hamming-weight-1 values
        std::vector<BinWord> preimages(input_length, 0);
        Integer n_found = 0;
        BinWord new_cstte = 0;
        for (BinWord x=1; x<(1 << input_length); x++)
        {
            BinWord y = cpp_linear_combination(image_vectors, x);
            if (cpp_hamming_weight(y) == 1)
            {
                unsigned int m = cpp_msb(y);
                if (preimages[m] != 0)
                    throw std::runtime_error("trying to invert non-invertible cpp_F2AffineMap");
                else
                {
                    preimages[m] = x;
                    n_found ++;
                }
            }
            if ((cstte != 0) and (y == cstte))
                new_cstte = x;
            if ((n_found == input_length) and ((cstte == 0) or (new_cstte != 0)))
                return cpp_F2AffineMap(preimages, input_length, input_length, new_cstte);
        }
        // this point can't be reached if the mapping is invertible
        throw std::runtime_error("trying to invert non-invertible cpp_F2AffineMap");
    }
    else
        // !TODO! implement inversion of cpp_F2AffineMapping in the general case
        throw std::runtime_error("inversion of F2AffineMap not implemented yet");
}


cpp_F2AffineMap cpp_F2AffineMap::transpose() const
{
    if (cstte != 0)
        throw std::runtime_error("trying to transpose non-linear cpp_F2AffineMap");
    else
    {
        std::vector<BinWord> result_imgs(output_length, 0);
        for(unsigned int i=0; i<output_length; i++)
        {
            BinWord mask = 1 << i;
            for(unsigned int j=0; j<input_length; j++)
                if (image_vectors[j] & mask)
                    result_imgs[i] ^= (1 << j) ;
        }
        return cpp_F2AffineMap(result_imgs, output_length, input_length);
    }
}


cpp_F2AffineMap operator+(const BinWord cst, const cpp_F2AffineMap & A) {
    return A + cst;
}

cpp_F2AffineMap identity_F2AffineMap(unsigned int n)
{
    std::vector<BinWord> imgs(n, 0);
    for(unsigned int i=0; i<n; i++)
        imgs[i] = ((BinWord)1) << i;
    return cpp_F2AffineMap(imgs, n, n);
}


cpp_F2AffineMap cpp_F2AffineMap_from_lut(Lut &s)
{
    std::vector<BinWord> images;
    for(unsigned int x=1; x<s.size(); x <<= 1)
        images.push_back(s[x] ^ s[0]);
    return cpp_F2AffineMap(images, s[0]);
}


cpp_F2AffineMap cpp_F2AffineMap_from_lut(cpp_S_box & s)
{
    std::vector<BinWord> images;
    for(unsigned int x=1; x<s.size(); x <<= 1)
        images.push_back(s[x] ^ s[0]);
    return cpp_F2AffineMap(images, s[0]);
}


// !TODO! A cpp_block_F2AffineMap that works to generate non-diagonal block matrices would be nice

cpp_F2AffineMap cpp_block_diagonal_F2AffineMap(
    const cpp_F2AffineMap &A,
    const cpp_F2AffineMap &B
    )
{
    std::vector<BinWord> images(B.image_vectors);
    for(auto v : A.image_vectors)
        images.push_back(v << B.get_output_length());
    return cpp_F2AffineMap(
        images,
        B.cstte | (A.cstte << B.get_output_length())
        );
}

// !TODO! Check with Léo that it is grosboutiste

/// @brief Computes the 2n by 2n affine map composed of the blocks A,B,C,D of the shape [[A,D][C,B]]. Assumes that all the blocks are squares of the same size 
/// @param A a n by n cpp_F2AffineMap
/// @param B a n by n cpp_F2AffineMap
/// @param C a n by n cpp_F2AffineMap
/// @param D a n by n cpp_F2AffineMap
/// @return a 2n by 2n cpp_F2AffineMap [[A,D][C,B]]
cpp_F2AffineMap cpp_F2AffineMap_from_blocks(
    const cpp_F2AffineMap &A,
    const cpp_F2AffineMap &B,
    const cpp_F2AffineMap &C,
    const cpp_F2AffineMap &D
    )
{

    if (A.get_input_length() != A.get_output_length()) 
        throw std::runtime_error("in cpp_EA_F2AffineMap: A cannot be invertible");
    else if (B.get_input_length() != B.get_output_length()) 
        throw std::runtime_error("in cpp_EA_F2AffineMap: B cannot be invertible");
    else if (A.get_input_length() < C.get_input_length())
        throw std::runtime_error("in cpp_EA_F2AffineMap: A and C have incompatible input length");
    else if (B.get_output_length() < C.get_output_length())
        throw std::runtime_error("in cpp_EA_F2AffineMap: B and C have incompatible output length");
    else if (A.get_input_length() < D.get_input_length())
        throw std::runtime_error("in cpp_EA_F2AffineMap: A and D have incompatible input length");
    else if (B.get_output_length() < D.get_output_length())
        throw std::runtime_error("in cpp_EA_F2AffineMap: B and D have incompatible output length");
    else
    {

        std::vector<BinWord> images;
        for(unsigned int i=0; i<B.get_input_length(); i++)
            images.push_back((D.image_vectors[i] << B.get_output_length()) |B.image_vectors[i]);
        for(unsigned int i=0; i<A.get_output_length(); i++)
            images.push_back((A.image_vectors[i] << B.get_output_length()) | C.image_vectors[i]);
        return cpp_F2AffineMap(images, cpp_oplus(B.cstte, C.cstte)|cpp_oplus(A.cstte, D.cstte) << B.get_output_length());

    }
}


