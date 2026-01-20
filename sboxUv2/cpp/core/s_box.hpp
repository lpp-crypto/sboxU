#ifndef _S_BOX_
#define _S_BOX_

#include "../common.hpp"
#include "prng.hpp"
#include "f2functions.hpp"

class cpp_Spectrum;

/** @class cpp_S_box

    This class stores the lookup table of a vectorial Boolean function and provides convenient methods to simplify the interaction with it. 
 */
class cpp_S_box
{
private:
    std::vector<BinWord> lut;
    Integer input_length;
    Integer output_length;

public:
    inline cpp_S_box() : lut(0), input_length(0), output_length(0) {} ;

    inline cpp_S_box(const cpp_S_box & s) :
        lut(s.get_lut()),
        input_length(s.get_input_length()),
        output_length(s.get_output_length())
    {} ;

    inline cpp_S_box(cpp_S_box * s) :
        lut(s->get_lut()),
        input_length(s->get_input_length()),
        output_length(s->get_output_length())
    {} ;

    cpp_S_box(std::vector<BinWord> _lut,
            Integer _input_length,
            Integer _output_length) ;

    cpp_S_box(std::vector<BinWord> _lut) ;

    cpp_S_box(Bytearray bytes) ;

    void destruct()
    {
        lut.clear();
        lut.shrink_to_fit();
    }
    
    ~cpp_S_box()
    {
        destruct();
    }

    /** Returns the lookup table of the binary sum (XOR) of this S-box and `s`, the sum being done for each possible input.
        
        @param s The other S-box.
        @return A new S-box containing the lookup of the XOR of the S-boxes.
     */
    cpp_S_box operator+ (const cpp_S_box &s) const;

    cpp_S_box operator* (const cpp_S_box &s) const;

    bool operator==(const cpp_S_box & other_s) const;

    std::string content_string_repr() const;

    Bytearray to_bytes() const;

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

    /** Returns a cpp_S_box implementing the compositional inverse of this one.

        If this S-box is not invertible, returns an empty
        S-box. Throwing an exception might appear cleaner, but in
        practice we couldn't convince cython to properly catch C++
        exceptions.
     */
    cpp_S_box inverse() const;

    cpp_S_box coordinate(BinWord i) const;
    
    cpp_S_box component(BinWord a) const;
    
    cpp_S_box derivative(BinWord delta) const;

};


cpp_S_box cpp_translation(const BinWord a, const Integer input_bit_length);

cpp_S_box cpp_empty_S_box();

Lut cpp_inverse(Lut & s);

bool cpp_is_permutation(Lut & s);

std::vector<BinWord> cpp_anf_component( const cpp_S_box & f);

cpp_Spectrum cpp_degree_spectrum(const cpp_S_box &f);


#endif
