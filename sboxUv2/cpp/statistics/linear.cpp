#include "linear.hpp"



// !SECTION! The LAT itself 

/* The algorithm (and its notations) come from:
 * 
 * Antoine Joux, "Algorithmic Cryptanalysis", Chapman & Hall/CRC,
 * 2009; Algorithm 9.3 (page 276).
 */
std::vector<Integer> cpp_walsh_transform(const cpp_S_box & f)
{
    // generating extended Boolean function
    std::vector<Integer> T(f.input_space_size(), 1);
    for (unsigned int x=0; x<f.input_space_size(); x++)
        if (f[x] == 1)
            T[x] = -1;
    // Walsh transform itself
    Integer sum = 0, diff = 0;
    for (unsigned int Sz=1; Sz<f.input_space_size(); Sz <<= 1)
    {
        for (unsigned int pos=0; pos<f.input_space_size(); pos += 2*Sz)
        {
            for (unsigned int j=0; j<Sz; j++)
            {
                unsigned int p = pos + j;
                sum  = T[p] + T[p + Sz];
                diff = T[p] - T[p + Sz];
                T[p] = sum;
                T[p + Sz] = diff;
            }
        }
    }
    return T;
}


std::vector< std::vector<Integer> > cpp_lat(const cpp_S_box & s)
{
    std::vector<std::vector<Integer> > table;
    table.reserve(s.output_space_size());
    // generating table
    for (unsigned int b=0; b<s.output_space_size(); b++)
        table.push_back(cpp_walsh_transform(s.component(b)));
    // transposing it
    std::vector<std::vector<Integer> > result;
    result.reserve(s.input_space_size());
    for (unsigned int a=0; a<s.input_space_size(); a++)
    {
        std::vector<Integer> row(s.output_space_size(), 0);
        for (unsigned int b=0; b<s.output_space_size(); b++)
            row[b] = table[b][a];
        result.push_back(row);
    }
    return result;
}


// !SECTION! Walsh spectrum


cpp_Spectrum cpp_walsh_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads)
{
    cpp_Spectrum count;
    int threads = threads_from_size(s.input_space_size());
#pragma omp parallel for reduction(aggregateSpectrum:count) num_threads(threads)
    for( unsigned int b = 1; b < s.input_space_size(); b++)
        count.incr_by_counting(cpp_walsh_transform(s.component(b)));
    return count;
}


cpp_Spectrum cpp_absolute_walsh_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads)
{
    return cpp_walsh_spectrum(s, n_threads).absolute();
}


cpp_S_box cpp_invert_lat(const std::vector< std::vector<Integer> > & l)
{
    std::vector<BinWord> result(l.size(), 0);
    for (unsigned int i=0; l[0].size() > (1 << i); i++)
    {
        Integer b = (1 << i), sum;
        for (unsigned int x=0; x<l.size(); x++)
        {
            sum = 0;
            for (unsigned int a=0; a<l.size(); a++)
            {
                if (cpp_scal_prod(a, x) == 0)
                    sum += l[a][b];
                else
                    sum -= l[a][b];
            }
            if (sum < 0)
                result[x] |= (1 << i) ;
        }
    }
    return cpp_S_box(result);
}
