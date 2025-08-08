#include "differential.hpp"



// !SECTION! The DDT itself 

std::vector<Integer> cpp_ddt_row(const cpp_S_box s, const BinWord delta)
{
    std::vector<Integer> result(s.output_space_size(), 0);
    for (unsigned int x=0; x<s.input_space_size(); x++)
        result[s[x^delta] ^ s[x]] ++ ;
    return result;
};



std::vector< std::vector<Integer> > cpp_ddt(const cpp_S_box s)
{
    std::vector< std::vector<Integer> > table ;
    table.reserve(s.input_space_size());
    for (unsigned int delta=0 ; delta<s.input_space_size(); delta++)
        table.push_back(cpp_ddt_row(s, delta));
    return table ;
}



// !SECTION! Differential spectrum 


void cpp_ddt_rows_count(
    cpp_Spectrum &result,
    const cpp_S_box s,
    const BinWord a_min,
    const BinWord a_max)
{
    for (unsigned int a=a_min; a<a_max; a++)
        result.incr_by_counting(cpp_ddt_row(s, a));
}



cpp_Spectrum cpp_differential_spectrum(
    const cpp_S_box s,
    const unsigned int n_threads)
{
    cpp_Spectrum count;
    if (n_threads == 1)
    {
        // small S-Box
        cpp_ddt_rows_count(std::ref(count),
                           s,
                           1,
                           s.input_space_size());
    }
    else
    {
        std::vector<std::thread> threads;
        std::vector<cpp_Spectrum> local_counts(n_threads);
        BinWord lower_bound = 1;
        for (unsigned int i=0; i<n_threads; i++)
        {
            // Will break on 32-bit arch if nthreads*s.size >= 1 << 32
            BinWord upper_bound = ((i+1)*s.input_space_size())/n_threads;
            threads.push_back(std::thread(cpp_ddt_rows_count,
                                          std::ref(local_counts[i]),
                                          s,
                                          lower_bound,
                                          upper_bound));
            lower_bound = upper_bound;

        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            count += local_counts[i];
        }
    }
    return count;
}



// !SECTION! Testing differential uniformity



bool cpp_is_ddt_row_max_smaller_than_2(
    const cpp_S_box s,
    const BinWord a)
{
    std::vector<uint64_t> row(s.input_space_size() >> 6,0);
    for (unsigned int x=0; x<s.input_space_size(); x++)
    {
        BinWord y = x^a;
        if (y < x)
            continue;
        BinWord d_out = s[y] ^ s[x];
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
    const cpp_S_box s,
    const BinWord a,
    const Integer u)
{
    std::vector<Integer> row(s.input_space_size(), 0);
    for (unsigned int x=0; x<s.input_space_size(); x++)
    {
        BinWord d_out = s[x^a] ^ s[x];
        row[d_out] ++ ;
        if (row[d_out] > u)
            return false;
    }
    return true;
}


bool cpp_is_ddt_row_max_smaller_than(
    const cpp_S_box s,
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


bool cpp_is_differential_uniformity_smaller_than_2(const cpp_S_box s)
{
    for (unsigned int a=1; a<s.input_space_size(); a++)
        if (cpp_is_ddt_row_max_smaller_than_2(s, a) == false)
            return false;
    return true;
}


bool cpp_is_differential_uniformity_smaller_than_u(
    const cpp_S_box s,
    const Integer u)
{
    for (unsigned int a=1; a<s.input_space_size(); a++)
        if (cpp_is_ddt_row_max_smaller_than_u(s, a, u) == false)
            return false;
    return true;
}


bool cpp_is_differential_uniformity_smaller_than(
    const cpp_S_box s,
    const Integer u)
{
    if(u == 2)
        return cpp_is_differential_uniformity_smaller_than_2(s);
    else
        return cpp_is_differential_uniformity_smaller_than_u(s,u);
}
