#include "./binLinearMap.hpp"



cpp_BinLinearMap cpp_BinLinearMap::operator*(const cpp_BinLinearMap l) const
{
    if (l.get_output_length() != input_length)
        throw std::runtime_error("mismatched input/output length when multiplying BinLinearMap:s");
    else
    {
        std::vector<BinWord> images(output_length, 0);
        for(unsigned int i=0; i<images.size(); i++)
            images[i] = cpp_linear_combination(image_vectors, l.image_vectors[i]);
        return cpp_BinLinearMap(images, input_length, output_length);
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
    // !TODO! implement inversion of BinLinearMapping 
    throw std::runtime_error("inversion of BinLinearMap not implemented yet");
}


cpp_BinLinearMap identity_BinLinearMap(unsigned int n)
{
    std::vector<BinWord> imgs(n, 0);
    for(unsigned int i=0; i<n; i++)
        imgs[i] = ((BinWord)1) << i;
    return cpp_BinLinearMap(imgs, n, n);
}
