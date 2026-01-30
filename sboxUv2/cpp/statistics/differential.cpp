#include "differential.hpp"

// !SECTION! The DDT itself 

std::vector<Integer> cpp_ddt_row(const cpp_S_box & s, const BinWord delta)
{
    std::vector<Integer> result(s.output_space_size(), 0);
    FOR_ENUMERATE_DIFFERENCE_COSETS(x,delta,s.input_space_size())
    {
        result[s[x] ^ s[x^delta]] += 2 ;
    }
    return result;
};





std::vector< std::vector<Integer>> cpp_ddt(const cpp_S_box & s)
{
    std::vector< std::vector<Integer>> table ;
    table.reserve(s.input_space_size());
    table.push_back(std::vector<Integer>(s.output_space_size(),0));
    table[0][0] = s.input_space_size();
    for (unsigned int delta=1 ; delta<s.input_space_size(); delta++)
        table.push_back(cpp_ddt_row(s, delta));
    return table ;
}

// !SECTION! The XDDT and its associated sets
// !SUBSECTION! The sets
std::vector< std::vector<BinWord>> cpp_xddt_row(const cpp_S_box & s, const BinWord delta)
{   std::vector< std::vector<BinWord>> result(s.output_space_size());
    BinWord y;
    for (unsigned int i=1; i< s.output_space_size(); i++)
        result[i].reserve(s.get_input_length()*2);
    FOR_ENUMERATE_DIFFERENCE_COSETS(x,delta,s.input_space_size()){
        y = s[x] ^ s[x^delta];
        result[y].push_back(x);
        result[y].push_back(x^delta);
    }
        
    return result;
}

std::vector<BinWord> cpp_xddt_entry(const cpp_S_box &s, const BinWord a, const BinWord b)
{   
    std::vector<BinWord> result;
    FOR_ENUMERATE_DIFFERENCE_COSETS(x,a,s.input_space_size()){
        if ((s[x]^s[x^a])==b){
            result.push_back(x);
            result.push_back(x^a);
        }
    }
    return result;
}

 Xtable cpp_xddt(const cpp_S_box & s){
    Xtable table; 
    table.reserve(s.input_space_size());

    std::vector<std::vector<BinWord>> first_row(s.input_space_size()); 
    for (BinWord x = 0; x < s.input_space_size(); x++) {
        first_row[0].push_back(x);
    }
    table.push_back(first_row);
    for (BinWord delta = 1; delta < s.input_space_size(); delta++) {
        table.push_back(cpp_xddt_row(s, delta));
}

return table;
}

std::vector<std::vector<BinWord>> cpp_yddt_row(const cpp_S_box & s, const BinWord delta)
{   std::vector< std::vector<BinWord>> result(s.output_space_size());
    BinWord y;
    for (unsigned int i=1; i< s.output_space_size(); i++)
        result[i].reserve(s.get_input_length()*2);
    FOR_ENUMERATE_DIFFERENCE_COSETS(x,delta,s.input_space_size()){
        y= s[x] ^ s[x^delta];
        result[y].push_back(s[x]);
        result[y].push_back(s[x^delta]);
    }
        
    return result;
}

 Xtable cpp_yddt(const cpp_S_box & s){
    Xtable table; 
    table.reserve(s.input_space_size());

    std::vector<std::vector<BinWord>> first_row(s.input_space_size()); 
    for (BinWord x = 0; x < s.input_space_size(); x++) {
        first_row[0].push_back(x);
    }
    table.push_back(first_row);
    for (BinWord delta = 1; delta < s.input_space_size(); delta++) {
        table.push_back(cpp_yddt_row(s, delta));
}

return table;
}


std::vector< std::vector<BinWord>> cpp_zddt_row(const cpp_S_box & s, const BinWord delta)
{   std::vector< std::vector<BinWord>> result(s.output_space_size());
    BinWord y;
    for (unsigned int i=1; i< s.output_space_size(); i++)
        result[i].reserve(s.get_input_length()*2);
    FOR_ENUMERATE_DIFFERENCE_COSETS(x,delta,s.input_space_size()){
        y = s[x] ^ s[x^delta];
        result[y].push_back((x|(s[x]<< s.get_input_length())));
        result[y].push_back((x^delta|(s[x^delta]<< s.get_input_length())));
    }
        
    return result;
}

 Xtable cpp_zddt(const cpp_S_box & s){
    Xtable table; 
    table.reserve(s.input_space_size());

    std::vector<std::vector<BinWord>> first_row(s.input_space_size()); 
    for (BinWord x = 0; x < s.input_space_size(); x++) {
        first_row[0].push_back(x|(s[x]<< s.get_input_length()));
    }
    table.push_back(first_row);
    for (BinWord delta = 1; delta < s.input_space_size(); delta++) {
        table.push_back(cpp_zddt_row(s, delta));
}

return table;
}




// !SECTION! Differential spectrum 


void cpp_ddt_rows_count(
    cpp_Spectrum & result,
    const cpp_S_box & s,
    const BinWord a_min,
    const BinWord a_max)
{
    for (unsigned int a=a_min; a<a_max; a++)
        result.incr_by_counting(cpp_ddt_row(s, a));
}


cpp_Spectrum cpp_differential_spectrum(const cpp_S_box & s,
    const unsigned int n_threads)
{
    cpp_Spectrum count;
    int threads = threads_from_size(s.input_space_size());
#pragma omp parallel for reduction(aggregateSpectrum:count) num_threads(threads)
    for( unsigned int a = 1; a < s.input_space_size(); a++)
        count.incr_by_counting(cpp_ddt_row(s,a));
    return count;
}


// !SECTION! Testing differential uniformity


bool cpp_is_ddt_row_max_smaller_than_2(
    const cpp_S_box & s,
    const BinWord a)
{
    std::vector<uint64_t> row((s.input_space_size() + 63) >> 6,0);

    FOR_ENUMERATE_DIFFERENCE_COSETS(x,a,s.input_space_size())
    {
        BinWord d_out = s[x] ^ s[x^a];
        // Assumes BinWord unsigned
        size_t index = d_out >> 6;
        uint64_t field = ((uint64_t) 1) << (d_out & 0x3F);
        if (row[index] & field)
            return false;
        row[index] ^= field;
    }
    return true;
}


bool cpp_is_ddt_row_max_smaller_than_u(
    const cpp_S_box & s,
    const BinWord a,
    const Integer u)
{
    std::vector<Integer> row(s.input_space_size(), 0);

    FOR_ENUMERATE_DIFFERENCE_COSETS(x,a,s.input_space_size())
    {
        BinWord d_out = s[x] ^ s[x^a];
        row[d_out] += 2 ;
        if (row[d_out] > u)
            return false;
    }
    return true;
}


bool cpp_is_ddt_row_max_smaller_than(
    const cpp_S_box & s,
    const BinWord a,
    const Integer u)
{
    if (a == 0)
        return s.input_space_size() <= u;
    else if (a >= s.input_space_size())
        // Should probably throw
        return true;
    else if (u == 2)
        return cpp_is_ddt_row_max_smaller_than_2(s,a);
    else
        return cpp_is_ddt_row_max_smaller_than_u(s,a,u);
}


bool cpp_is_differential_uniformity_smaller_than_2(const cpp_S_box & s)
{
    for (unsigned int a=1; a<s.input_space_size(); a++)
        if (cpp_is_ddt_row_max_smaller_than_2(s, a) == false)
            return false;
    return true;
}


bool cpp_is_differential_uniformity_smaller_than_u(
    const cpp_S_box & s,
    const Integer u)
{
    for (unsigned int a=1; a<s.input_space_size(); a++)
        if (cpp_is_ddt_row_max_smaller_than_u(s, a, u) == false)
            return false;
    return true;
}


bool cpp_is_differential_uniformity_smaller_than(
    const cpp_S_box & s,
    const Integer u)
{
    if(u == 2)
        return cpp_is_differential_uniformity_smaller_than_2(s);
    else
        return cpp_is_differential_uniformity_smaller_than_u(s,u);
}
