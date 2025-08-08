#include "../f2functions.hpp"
#include "linear.hpp"



// !SECTION! The LAT itself 

/* The algorithm (and its notations) come from:
 * 
 * Antoine Joux, "Algorithmic Cryptanalysis", Chapman & Hall/CRC,
 * 2009; Algorithm 9.3 (page 276).
 */
std::vector<Integer> cpp_walsh_transform(const cpp_S_box f)
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


std::vector< std::vector<Integer> > cpp_lat(const cpp_S_box s)
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


void walsh_spectrum_cols_count(
    cpp_Spectrum &result,
    const cpp_S_box s,
    const BinWord b_min,
    const BinWord b_max)
{
    for (unsigned int b=b_min; b<b_max; b++)
        result.incr_by_counting(cpp_walsh_transform(s.component(b))) ;
}


cpp_Spectrum cpp_walsh_spectrum(
    const cpp_S_box s,
    const unsigned int n_threads)
{
    cpp_Spectrum count;
    if (n_threads == 1)         // single thread
    {
        walsh_spectrum_cols_count(std::ref(count), s, 1, s.size());
    }
    else                        // strictly more than 1 thread
    {
        std::vector<std::thread> threads;
        std::vector<cpp_Spectrum> local_counts(n_threads, cpp_Spectrum());
        BinWord lower_bound = 1;
        for (unsigned int i=0; i<n_threads; i++)
        {
            // Will break on 32-bit arch if nthreads*s.size >= 1 << 32
            BinWord upper_bound = ((i+1)*s.output_space_size())/n_threads;
            threads.push_back(std::thread(walsh_spectrum_cols_count,
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


cpp_S_box cpp_invert_lat(const std::vector< std::vector<Integer> > l)
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


// // !SUBSECTION! Internal routines 




// void lat_zeroes_in_columns(
//     std::vector<BinWord> &result,
//     const Sbox s,
//     const BinWord b_min,
//     const BinWord b_max,
//     const unsigned int n)
// {
//     result.reserve(result.size() + b_max - b_min) ;
//     for (unsigned int b=b_min; b<b_max; b++)
//     {
//         // computing one coordinate
//         Sbox f(component_cpp(b, s)) ;
//         // Walsh transform
//         std::vector<Integer> w(walsh_spectrum_coord(f));
//         // getting zeroes
//         for (unsigned int a=0; a<w.size(); a++)
//             if (w[a] == 0)
//                 result.push_back((a << n) | b) ;
//     }
// }


// void do_lat_columns_contain_zero(
//     std::vector<bool> &result,
//     const Sbox s,
//     const BinWord b_min,
//     const BinWord b_max)
// {
//     result.reserve(result.size() + b_max - b_min) ;
//     for (unsigned int b=b_min; b<b_max; b++)
//     {
//         // computing one coordinate
//         Sbox f(component_cpp(b, s)) ;
//         // Walsh transform
//         std::vector<Integer> w(walsh_spectrum_coord(f));
//         bool contains_zero = false;
//         for (auto & c : w)
//             if (c == 0)
//             {
//                 contains_zero = true;
//                 break;
//             }
//         result.push_back(contains_zero) ;
//     }
// }

// // !SUBSECTION! High level parallelized functions 

// std::vector<BinWord> lat_zeroes_cpp(
//     const Sbox s,
//     const unsigned int n,
//     const unsigned int n_threads)
// {
//     check_length_cpp(s);
//     std::vector<BinWord> zeroes;
//     if (n_threads == 1)
//     {
//         // small S-Box
//         lat_zeroes_in_columns(std::ref(zeroes), s, 0, s.size(), n);
//     }
//     else
//     {
//         std::vector<std::thread> threads;
//         std::vector<std::vector<BinWord> > local_zeroes(n_threads);
//         unsigned int slice_size = s.size()/n_threads;
//         for (unsigned int i=0; i<n_threads; i++)
//         {
//             unsigned int
//                 lower_bound = i*slice_size,
//                 upper_bound = (i+1)*slice_size;
//             if (upper_bound > s.size())
//                 upper_bound = s.size();
//             threads.push_back(std::thread(lat_zeroes_in_columns,
//                                           std::ref(local_zeroes[i]),
//                                           s,
//                                           lower_bound,
//                                           upper_bound,
//                                           n));

//         }
//         for (unsigned int i=0; i<n_threads; i++)
//         {
//             threads[i].join();
//             zeroes.insert(zeroes.end(),
//                           local_zeroes[i].begin(),
//                           local_zeroes[i].end());
//         }
//     }
//     return zeroes ;
// }


// std::vector<BinWord> projected_lat_zeroes_cpp(
//     const Sbox s,
//     const unsigned int n_threads)
// {
//     check_length_cpp(s);
//     std::vector<bool> projection;
//     if (n_threads == 1)
//     {
//         // small S-Box
//         do_lat_columns_contain_zero(std::ref(projection), s, 0, s.size());
//     }
//     else
//     {
//         std::vector<std::thread> threads;
//         std::vector<std::vector<bool> > local_projections(n_threads);
//         unsigned int slice_size = s.size()/n_threads;
//         for (unsigned int i=0; i<n_threads; i++)
//         {
//             unsigned int
//                 lower_bound = i*slice_size,
//                 upper_bound = (i+1)*slice_size;
//             if (upper_bound > s.size())
//                 upper_bound = s.size();
//             threads.push_back(std::thread(do_lat_columns_contain_zero,
//                                           std::ref(local_projections[i]),
//                                           s,
//                                           lower_bound,
//                                           upper_bound));

//         }
//         for (unsigned int i=0; i<n_threads; i++)
//         {
//             threads[i].join();
//             projection.insert(projection.end(),
//                           local_projections[i].begin(),
//                           local_projections[i].end());
//         }
//     }
//     std::vector<BinWord> projected_zeroes;
//     for (unsigned int i=1; i<projection.size(); i++) // we purposefully leave out the zero
//         if (projection[i])
//             projected_zeroes.push_back(i) ;
//     return projected_zeroes ;
// }
