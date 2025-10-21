#ifndef _CORE_BINLINEARMAP_
#define _CORE_BINLINEARMAP_

#include "../sboxU.hpp"
#include "./f2functions.hpp"
#include "./s_box.hpp"

class cpp_BinLinearMap
{
private:
    Integer input_length;
    Integer output_length;

public:
    std::vector<BinWord> image_vectors;
    
    inline cpp_BinLinearMap() :
        image_vectors(std::vector<BinWord>(0)),
        input_length(0),
        output_length(0)
    {};
    
    
    inline cpp_BinLinearMap(const std::vector<BinWord> & _image_vectors) :
        image_vectors(_image_vectors.cbegin(), _image_vectors.cend()),
        input_length(_image_vectors.size()),
        output_length(0)
    {
        for(auto &v : image_vectors)
        {
            BinWord m = cpp_msb(v);
            if (m > output_length)
                output_length = m;
        }
        output_length ++;       // cpp_msb starts at 0
    };

    
    inline cpp_BinLinearMap(
        const std::vector<BinWord> & _image_vectors,
        const Integer _input_length,
        const Integer _output_length
        ) :
        image_vectors(_image_vectors),
        input_length(_input_length),
        output_length(_output_length)
    {};

    
    cpp_BinLinearMap(const cpp_S_box & lut) ;

    
    inline Integer get_input_length() const
    {
        return input_length;
    };
    
    
    inline Integer get_output_length() const
    {
        return output_length;
    };


        
    inline BinWord operator()(const BinWord x) const
    {
        return cpp_linear_combination(image_vectors, x);
    };

    
    cpp_BinLinearMap operator*(const cpp_BinLinearMap l) const;
    
    cpp_BinLinearMap operator+(const cpp_BinLinearMap l) const;

    cpp_BinLinearMap inverse() const;

    
    cpp_BinLinearMap transpose() const;


    inline BinWord rank() const
    {
        return cpp_rank_of_vector_set(image_vectors);
    };
    
    
    inline cpp_S_box get_cpp_S_box() const
    {
        std::vector<BinWord> lut(1 << input_length, 0);
        for(BinWord x=1; x<lut.size(); x++)
            lut[x] = cpp_linear_combination(image_vectors, x);
        return cpp_S_box(lut, input_length, output_length);
    };

    
    inline std::vector<BinWord> get_image_vectors() const
    {
        std::vector<BinWord> result(image_vectors);
        return result;
    };
    
};

cpp_BinLinearMap identity_BinLinearMap(unsigned int n);

cpp_BinLinearMap cpp_BinLinearMap_from_lut(Lut & s);

cpp_BinLinearMap cpp_BinLinearMap_from_lut(cpp_S_box & s);

cpp_BinLinearMap cpp_block_diagonal_BinLinearMap(
    const cpp_BinLinearMap &A,
    const cpp_BinLinearMap &B
    );

    
#endif
