#include "sboxU.hpp"



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

// !SECTION! Operator overloading 

BinWord cpp_S_box::operator[] (const BinWord x) const
{
    return lut[x];
}


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


// !SECTION! Basic access methods


std::vector<BinWord> cpp_S_box::get_lut() const
{
    return lut;
}


std::string cpp_S_box::content_string_repr() const
{
    std::stringstream result;
    result << "[";
    for (unsigned int x=0; x<size()-1; x++)
        result << std::hex << lut[x] << ",";
    result << std::hex << lut.back() << "]";
    return result.str();    
}


// !SUBSECTION! Sizes


Integer cpp_S_box::size() const
{
    return lut.size();
}

Integer cpp_S_box::get_input_length() const
{
    return input_length;
}


Integer cpp_S_box::input_space_size() const
{
    return (1 << input_length);
}

Integer cpp_S_box::get_output_length() const
{
    return output_length;
}


Integer cpp_S_box::output_space_size() const
{
    return (1 << output_length);
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
    for(unsigned int x=0; x<result_lut.size(); x++)
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
