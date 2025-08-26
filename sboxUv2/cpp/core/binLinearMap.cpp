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
    if (input_length <= 8)      // up to n=8, n^2.8 is greater than 2^n
    {
        // bruteforcing inputs and sotring those mapped to
        // hamming-weigh-1 values
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
                return cpp_BinLinearMap(preimages);
        }
        // this point can't be reached if the mapping is invertible
        throw std::runtime_error("trying to invert non-invertible cpp_BinLinearMap");
    }
    else
        // !TODO! implement inversion of cpp_BinLinearMapping in the general case
        throw std::runtime_error("inversion of BinLinearMap not implemented yet");
}


cpp_BinLinearMap identity_BinLinearMap(unsigned int n)
{
    std::vector<BinWord> imgs(n, 0);
    for(unsigned int i=0; i<n; i++)
        imgs[i] = ((BinWord)1) << i;
    return cpp_BinLinearMap(imgs, n, n);
}
