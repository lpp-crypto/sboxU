#include "./binLinearMap.hpp"


cpp_BinLinearMap::cpp_BinLinearMap(const cpp_S_box & lut) :
    image_vectors(0),
    input_length(0),
    output_length(0) 
{
    for(unsigned int pos=1; pos<lut.size(); pos<<=1)
    {
        BinWord
            y = lut[pos],
            m = cpp_msb(y);
        image_vectors.push_back(y);
        input_length += 1;
        if (m > output_length)
            output_length = m;
    }
    output_length ++;
}


cpp_BinLinearMap cpp_BinLinearMap::operator*(const cpp_BinLinearMap l) const
{
    if (l.get_output_length() != input_length)
        throw std::runtime_error("mismatched input/output length when multiplying BinLinearMap:s");
    else
    {
        std::vector<BinWord> images(l.get_input_length(), 0);
        for(unsigned int i=0; i<l.get_input_length(); i++)
            images[i] = cpp_linear_combination(image_vectors, l.image_vectors[i]);
        return cpp_BinLinearMap(images, l.get_input_length(), output_length);
    }
};


cpp_BinLinearMap cpp_BinLinearMap::operator+(const cpp_BinLinearMap l) const
{
    if (l.get_input_length() != input_length)
        throw std::runtime_error("mismatched input length when adding BinLinearMap:s");
    else
    {
        std::vector<BinWord> images(image_vectors.begin(),
                                    image_vectors.end());
        for(unsigned int i=0; i<images.size(); i++)
            images[i] ^= l.image_vectors[i];
        return cpp_BinLinearMap(images, input_length, output_length);
    }
};



cpp_BinLinearMap cpp_BinLinearMap::inverse() const
{
    if (input_length <= 20)      // up to n=8, n^2.8 is greater than 2^n
                                 // but for now, we don't have anything else
    {
        // bruteforcing inputs and sorting those mapped to
        // hamming-weight-1 values
        std::vector<BinWord> preimages(input_length, 0);
        Integer n_found = 0;
        for (BinWord x=1; x<(1 << input_length); x++)
        {
            BinWord y = cpp_linear_combination(image_vectors, x);
            if (cpp_hamming_weight(y) == 1)
            {
                unsigned int m = cpp_msb(y);
                if (preimages[m] != 0)
                    throw std::runtime_error("trying to invert non-invertible cpp_BinLinearMap");
                else
                {
                    preimages[m] = x;
                    n_found ++;
                }
            }
            if (n_found == input_length)
                return cpp_BinLinearMap(preimages, input_length, input_length);
        }
        // this point can't be reached if the mapping is invertible
        throw std::runtime_error("trying to invert non-invertible cpp_BinLinearMap");
    }
    else
        // !TODO! implement inversion of cpp_BinLinearMapping in the general case
        throw std::runtime_error("inversion of BinLinearMap not implemented yet");
}


cpp_BinLinearMap cpp_BinLinearMap::transpose() const
{
    std::vector<BinWord> result_imgs(output_length, 0);
    for(unsigned int i=0; i<output_length; i++)
    {
        BinWord mask = 1 << i;
        for(unsigned int j=0; j<input_length; j++)
            if (image_vectors[j] & mask)
                result_imgs[i] ^= (1 << j) ;
    }
    return cpp_BinLinearMap(result_imgs, output_length, input_length);
}


cpp_BinLinearMap identity_BinLinearMap(unsigned int n)
{
    std::vector<BinWord> imgs(n, 0);
    for(unsigned int i=0; i<n; i++)
        imgs[i] = ((BinWord)1) << i;
    return cpp_BinLinearMap(imgs, n, n);
}


cpp_BinLinearMap cpp_BinLinearMap_from_lut(Lut & s)
{
    if (s[0] != 0)
        throw std::runtime_error("Input of cpp_BinLinearMap_from_lut must map 0 to 0!");
    else
    {
        std::vector<BinWord> images;
        for(unsigned int x=1; x<s.size(); x <<= 1)
            images.push_back(s[x]);
        return cpp_BinLinearMap(images);
    }
}


cpp_BinLinearMap cpp_BinLinearMap_from_lut(cpp_S_box & s)
{
    if (s[0] != 0)
        throw std::runtime_error("Input of cpp_BinLinearMap_from_lut must map 0 to 0!");
    else
    {
        std::vector<BinWord> images;
        for(unsigned int x=1; x<s.size(); x <<= 1)
            images.push_back(s[x]);
        return cpp_BinLinearMap(images);
    }
}

// !TODO! A cpp_block_BinLinearMap that works to generate non-diagonal block matrices would be nice

cpp_BinLinearMap cpp_block_diagonal_BinLinearMap(
    const cpp_BinLinearMap &A,
    const cpp_BinLinearMap &B
    )
{
    std::vector<BinWord> images(B.image_vectors);
    for(auto v : A.image_vectors)
        images.push_back(v << B.get_output_length());
    return cpp_BinLinearMap(images);
}
