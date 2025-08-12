#ifndef _S_BOX_
#define _S_BOX_

#include "sboxU.hpp"
#include "f2functions.hpp"


class cpp_S_box
{
private:
    std::vector<BinWord> lut;
    Integer input_length;
    Integer output_length;

public:
    inline cpp_S_box() : lut(0), input_length(0), output_length(0) {} ;
    
    cpp_S_box(std::vector<BinWord> _lut,
            Integer _input_length,
            Integer _output_length) ;

    cpp_S_box(std::vector<BinWord> _lut) ;


    cpp_S_box operator+ (const cpp_S_box &s) const;

    cpp_S_box operator* (const cpp_S_box &s) const;

    std::string content_string_repr() const;

    inline BinWord operator[] (const BinWord x) const
    {
        return lut[x];
    };

    inline std::vector<BinWord> get_lut() const
    {
        return lut;
    };
    

    inline Integer size() const
    {
        return lut.size();
    };

    inline Integer get_input_length() const
    {
        return input_length;
    };

    inline Integer input_space_size() const
    {
        return (1 << input_length);
    };

    inline Integer get_output_length() const
    {
        return output_length;
    };
    
    Integer output_space_size() const
    {
        return (1 << output_length);
    };
    

    bool is_invertible() const;

    cpp_S_box inverse() const;

    cpp_S_box coordinate(BinWord i) const;
    
    cpp_S_box component(BinWord a) const;
    
    cpp_S_box derivative(BinWord delta) const;

    
    // !TODO! implement:
    // ! - ANF (but not in a BooleanFunction way)

    // !TODO! documentation:
    // ! - endianess
    // ! - swap_matrix and tu_decomposition are not coherent
};


cpp_S_box cpp_translation(const BinWord a, const Integer input_bit_length);

cpp_S_box cpp_empty_S_box();

#endif
