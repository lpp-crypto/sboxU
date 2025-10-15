#include "s_box.hpp"

// !SECTION! Constructors 



cpp_S_box::cpp_S_box(std::vector<BinWord> _lut,
                 Integer _input_length,
                 Integer _output_length) :
    lut(_lut),
    input_length(_input_length),
    output_length(_output_length)
{}



cpp_S_box::cpp_S_box(std::vector<BinWord> _lut) :
    lut(_lut),
    input_length(1),
    output_length(1)
{
    // setting output_length
    for (unsigned int x=0; x<lut.size(); x++)
    {
        BinWord m = cpp_msb(lut[x])+1;
        if (m > output_length)
            output_length = m ;
    }
    // setting input_length
    input_length = cpp_msb(lut.size());
}


cpp_S_box::cpp_S_box(Bytearray bytes) 
{
    output_length = bytes[0];
    if (output_length <= 4)
    {
        // case where we store two outputs per byte
        lut.assign((bytes.size()-1) * 2, 0);
        for(BinWord x=0; x<lut.size(); x+=2)
        {
            BinWord y = bytes[(x >> 1) + 1];
            lut[x]   = y & 0xF;
            lut[x+1] = y >> 4;
        }
    }
    else
    {
        // case where we store outputs over whole bytes
        unsigned int output_byte_length = output_length / 8;
        if ((output_length % 8) != 0)
            output_byte_length ++ ;
        lut.assign((bytes.size()-1) / output_byte_length, 0);
        for(BinWord x=0; x<lut.size(); x++)
            for(unsigned int i=0; i<output_byte_length; i++)
                lut[x] ^= bytes[x*output_byte_length + i + 1] << (i*8);
    }
    input_length = cpp_msb(lut.size());
}


// !SECTION! Operator overloading 

cpp_S_box cpp_S_box::operator+ (const cpp_S_box &s) const
{
    if (s.get_input_length() != input_length)
        throw std::runtime_error("Trying to add S_boxes of different sizes");
    else
    {
        Integer new_output_length = output_length;
        if (s.get_output_length() > new_output_length)
            new_output_length = s.get_output_length();
        std::vector<BinWord> new_lut(lut.begin(), lut.end());
        for (unsigned int x=0; x<lut.size(); x++)
            new_lut[x] ^= s[x];
        return cpp_S_box(new_lut, input_length, new_output_length);
    }
}


cpp_S_box cpp_S_box::operator* (const cpp_S_box &s) const
{
    if (input_length != s.get_output_length())
        throw std::runtime_error("Trying to compose S_boxes of incompatible sizes");
    else
    {
        std::vector<BinWord> new_lut(lut.size(), 0);
        for (unsigned int x=0; x<s.input_space_size(); x++)
            new_lut[x] = lut[s[x]];
        return cpp_S_box(new_lut, s.get_input_length(), output_length);
    }
}


bool cpp_S_box::operator==(const cpp_S_box & other_s) const
{
    if ((input_length != other_s.get_input_length())
        or
        (output_length != other_s.get_output_length()))
        return false;
    for(BinWord x=0; x<input_space_size(); x++)
        if (lut[x] != other_s[x])
            return false;
    return true;
}


// !SECTION! Basic access methods



std::string cpp_S_box::content_string_repr() const
{
    if (lut.size() > 0)
    {
        std::stringstream result;
        result << "[";
        for (unsigned int x=0; x<size()-1; x++)
            result << std::dec << lut[x] << ",";
        result << std::dec << lut.back() << "]";
        return result.str();    
    }
    else
        return "[]";
}


Bytearray cpp_S_box::to_bytes() const
{
    if (output_length <= 4)
    {
        // case where we store two outputs per byte
        Bytearray result(lut.size() / 2 + 1, 0);
        result[0] = output_length;
        for(BinWord x=0; x<lut.size(); x+=2)
        {
            BinWord x_prime = (x >> 1) + 1;
            result[x_prime] = lut[x] ;
            result[x_prime] ^= lut[x+1] << 4;
        }
        return result;
    }
    else
    {
        // case where each output is stored over whole bytes
        unsigned int output_byte_length = output_length / 8;
        if ((output_length % 8) != 0)
            output_byte_length ++ ;
        Bytearray result(lut.size() * output_byte_length + 1, 0);
        result[0] = output_length;
        for(BinWord x=0; x<lut.size(); x++)
            for(BinWord i=0; i<output_byte_length; i++)
                result[x*output_byte_length + i + 1] = (lut[x] >> (i*8)) & 0xFF;
        return result;
    }
}


// !SECTION! More sophisticated operations


// !SUBSECTION! Inversion 


bool cpp_S_box::is_invertible() const
{
    if (input_length != output_length)
        return false;
    else
    {
        std::vector<BinWord> counters(lut.size(), 0);
        for(unsigned int x=0; x<input_space_size(); x++)
            if ((counters[lut[x]]++) > 1)
                return false;
        return true;
    }
}

cpp_S_box cpp_S_box::inverse() const
{
    if (input_length != output_length)
        throw std::runtime_error("The S_box is not invertible");
    else
    {
        BinWord forbidden_value = 1 << (output_length + 1);
        std::vector<BinWord> inverse_lut(input_space_size(), forbidden_value);
        for(unsigned int x=0; x<input_space_size(); x++)
        {
            if (inverse_lut[lut[x]] != forbidden_value)
                throw std::runtime_error("The S_box is not invertible");
            else
                inverse_lut[lut[x]] = x;
        }
        return cpp_S_box(inverse_lut, input_length, input_length);
    }
}


// !SUBSECTION! Basic math operations


cpp_S_box cpp_S_box::coordinate(BinWord i) const
{
    std::vector<uint64_t> result_lut(input_space_size(), 0);
    for(unsigned int x=0; x<result_lut.size(); x++)
        result_lut[x] = (lut[x] >> i) & 1;
    return cpp_S_box(result_lut, input_length, 1);
}


cpp_S_box cpp_S_box::component(BinWord a) const
{
    std::vector<uint64_t> result_lut(input_space_size(), 0);
    for(unsigned int x=0; x<input_space_size(); x++)
        result_lut[x] = cpp_scal_prod(lut[x], a);
    return cpp_S_box(result_lut, input_length, 1);
}


cpp_S_box cpp_S_box::derivative(BinWord delta) const
{
    std::vector<uint64_t> result_lut(lut.begin(), lut.end());
    for(unsigned int x=0; x<result_lut.size(); x++)
        result_lut[x] ^= lut[x ^ delta];
    return cpp_S_box(result_lut, input_length, output_length);
}



// !SECTION! Other generating functions

cpp_S_box cpp_translation(const BinWord a, const Integer input_length)
{
    BinWord input_space_size = ((BinWord)1) << input_length;
    std::vector<BinWord> lut(input_space_size, 0);
    for(BinWord x=0; x<input_space_size; x++)
        lut[x] = x ^ a;
    return cpp_S_box(lut);
}

cpp_S_box cpp_empty_S_box()
{
    return cpp_S_box(std::vector<BinWord>(0), 0, 0);
}

// !SECTION! Basic helper functions

Lut cpp_inverse(Lut & s)
{
    Lut result(s.size(), 0);
    for (unsigned int x=0; x<s.size(); x++)
        result[s[x]] = x;
    return result;
}


bool cpp_is_permutation(Lut & s)
{
    std::vector<BinWord> counters(s.size(), 0);
    for(unsigned int x=0; x<s.size(); x++)
        if ((counters[s[x]]++) > 1)
            return false;
    return true;
}


