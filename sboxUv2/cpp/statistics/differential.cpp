#include "../sboxU.hpp"




std::vector<Integer> cpp_ddt_row(const cpp_S_box s, const BinWord a)
{
    std::vector<Integer> result(s.size(), 0);
    for (unsigned int x=0; x< s.size(); x++)
        result[s[x^a] ^ s[x]] ++ ;
    return result;
};


void cpp_ddt_rows_count(
    cpp_Spectrum &result,
    const cpp_S_box s,
    const BinWord a_min,
    const BinWord a_max)
{
    for (unsigned int a=a_min; a<a_max; a++)
        result.incr_by_counting(cpp_ddt_row(s, a));
}



std::vector< std::vector<Integer> > cpp_ddt(const cpp_S_box s)
{
    std::vector< std::vector<Integer> > table ;
    table.reserve(s.input_space_size());
    for (unsigned int i = 0 ; i < s.output_space_size(); i++)
        table.push_back(cpp_ddt_row(s, i));
    return table ;
}


cpp_Spectrum cpp_differential_spectrum(
    const cpp_S_box s,
    const unsigned int n_threads)
{
    cpp_Spectrum count;
    if (n_threads == 1)
    {
        // small S-Box
        cpp_ddt_rows_count(std::ref(count), s, 1, s.size());
    }
    else
    {
        std::vector<std::thread> threads;
        std::vector<cpp_Spectrum> local_counts(n_threads);
        BinWord lower_bound = 1;
        for (unsigned int i=0; i<n_threads; i++)
        {
            // Will break on 32-bit arch is nthreads*s.size >= 1 << 32
            BinWord upper_bound = ((i+1)*s.size())/n_threads;
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


