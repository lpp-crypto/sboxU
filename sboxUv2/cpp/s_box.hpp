#ifndef _SBOX_
#define _SBOX_


#include "sboxU.hpp"


class cpp_SBox
{
private:
    std::vector<BinWord> lut;
    Integer input_length;
    Integer output_length;

public:
    inline cpp_SBox() : lut(0), input_length(0), output_length(0) {} ;
    
    cpp_SBox(std::vector<BinWord> _lut,
            Integer _input_length,
            Integer _output_length) ;

    cpp_SBox(std::vector<BinWord> _lut) ;

    BinWord operator[] (const BinWord x) const;

    cpp_SBox operator+ (const cpp_SBox &s) const;

    cpp_SBox operator* (const cpp_SBox &s) const;

    std::string content_string_repr() const;

    std::vector<BinWord> get_lut() const;

    Integer size() const;

    Integer get_input_length() const;

    Integer input_space_size() const;

    Integer get_output_length() const;
    
    Integer output_space_size() const;


    bool is_invertible() const;

    cpp_SBox inverse() const;
    
    // !TODO! implement:
    // ! - .derivative(a)
    // ! - .component(a)
    // ! - ANF (but not in a BooleanFunction way)
    // ! - Component, Component of ANF


    // !TODO! documentation:
    // ! - endianess
    // ! - swap_matrix and tu_decomposition are not coherent
};


cpp_SBox cpp_translation(const BinWord a, const Integer input_bit_length);


#endif
