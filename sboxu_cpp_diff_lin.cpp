/* Time-stamp: <2017-10-06 19:03:02 lperrin>
 *
 * LICENSE
 */ 

#include "sboxu_cpp_diff_lin.hpp"
using namespace boost::python;


// !SECTION! Differential properties

// !SUBSECTION! Internal routines 

std::vector<Integer> ddt_row(const Sbox s, const uint32_t a)
{
    // I tried implementing this function using list
    // instead of vectors. However, it turned out to be about 2-3
    // times slower.
    std::vector<Integer> result(s.size(), 0);
    for (unsigned int x=0; x<s.size(); x++)
        result[s[x^a] ^ s[x]] ++ ;
    return result;
}


void ddt_rows_count(
    std::map<Integer, Integer> &result,
    const Sbox s,
    const uint32_t a_min,
    const uint32_t a_max)
{
    for (unsigned int a=a_min; a<a_max; a++)
    {
        std::vector<Integer> row(ddt_row(s, a));
        for (auto &v : row)
            result[v] ++;
    }
}


// !SUBSECTION! Python-facing functions 

list ddt(const list& l)
{
    list result;    
    Sbox s(lst_2_vec_int(l)) ;
    check_length(s);
    for (unsigned int a=0; a<s.size(); a++)
        result.append(vec_2_lst_int(ddt_row(s, a)));

    return result;
}


dict differential_spectrum_fast(const list& l, const unsigned int n_threads)
{
    Sbox s(lst_2_vec_int(l)) ;
    check_length(s);
    std::map<Integer, Integer> count;
    if (n_threads == 1)
    {
        // small S-Box
        ddt_rows_count(std::ref(count), s, 1, s.size());
    }
    else
    {
        std::vector<std::thread> threads;
        std::vector<std::map<Integer, Integer> > local_counts(n_threads);
        unsigned int slice_size = s.size()/n_threads;
        for (unsigned int i=0; i<n_threads; i++)
        {
            unsigned int
                lower_bound = i*slice_size,
                upper_bound = (i+1)*slice_size;
            if (lower_bound == 0)
                lower_bound = 1;
            if (upper_bound > s.size())
                upper_bound = s.size();
            threads.push_back(std::thread(ddt_rows_count,
                                          std::ref(local_counts[i]),
                                          s,
                                          lower_bound,
                                          upper_bound));

        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            for (auto &entry : local_counts[i])
                count[entry.first] += entry.second ;
        }
    }
    dict result;
    for (auto &entry : count)
        result[entry.first] = entry.second ;
    return result;
}


// !SECTION! Linear properties


// !SUBSECTION! Internal routines 

/* f is assumed to be the truth-table of a Boolean function. Its size
 * must thus be a power of two. At the end of the function, the output
 * T is the Walsh spectrum of f. If n is the size of the input of f,
 * then:
 *
 * T(a) = \\sum_{x \\in \\{ 0,1 \\}^n}} (-1)^{<a,x> + f(x)} .
 *
 * The algorithm (and its notations) come from:
 * 
 * Antoine Joux, "Algorithmic Cryptanalysis", Chapman & Hall/CRC,
 * 2009; Algorithm 9.3 (page 276).
 */
std::vector<Integer> walsh_spectrum_coord(const Sbox f)
{
    // generating extended Boolean function
    std::vector<Integer> T(f.size(), 1);
    for (unsigned int x=0; x<f.size(); x++)
        if (f[x] == 1)
            T[x] = -1;
    // Walsh transform itself
    Integer sum = 0, diff = 0;
    for (unsigned int Sz=1; Sz<f.size(); Sz <<= 1)
    {
        for (unsigned int pos=0; pos<f.size(); pos += 2*Sz)
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


void walsh_spectrum_cols_count(
    std::map<Integer, Integer> &result,
    const Sbox s,
    const uint32_t b_min,
    const uint32_t b_max)
{
    for (unsigned int b=b_min; b<b_max; b++)
    {
        // computing one coordinate
        Sbox f(s.size(), 0);
        for (unsigned int x=0; x<f.size(); x++)
            f[x] = scal_prod(b, s[x]);
        // Walsh transform
        std::vector<Integer> w(walsh_spectrum_coord(f));
        // counting
        for (auto &v : w)
            result[v] ++;
    }
}


// !SUBSECTION! Python-facing functions 

list lat(const list& l)
{
    Sbox s(lst_2_vec_int(l)) ;
    check_length(s);
    std::vector<std::vector<Integer> > table(
        s.size(),
        std::vector<Integer>(s.size(), 0));
    // generating table
    for (unsigned int b=0; b<s.size(); b++)
    {
        // computing one coordinate
        Sbox f(s.size(), 0);
        for (unsigned int x=0; x<f.size(); x++)
            f[x] = scal_prod(b, s[x]);
        // Walsh transform
        std::vector<Integer> w(walsh_spectrum_coord(f));
        table[b].assign(w.begin(), w.end()); 
    }
    // transposing
    list result;
    for (unsigned int a=0; a<s.size(); a++)
    {
        std::vector<Integer> row(s.size(), 0);
        for (unsigned int b=0; b<s.size(); b++)
            row[b] = table[b][a];
        result.append<list>(vec_2_lst_int(row));
    }
    return result;
}


dict walsh_spectrum_fast(const list& l, const unsigned int n_threads)
{
    Sbox s(lst_2_vec_int(l)) ;
    check_length(s);
    std::map<Integer, Integer> count;
    if (n_threads == 1)
    {
        // small S-Box
        walsh_spectrum_cols_count(std::ref(count), s, 1, s.size());
    }
    else
    {
        std::vector<std::thread> threads;
        std::vector<std::map<Integer, Integer> > local_counts(n_threads);
        unsigned int slice_size = s.size()/n_threads;
        for (unsigned int i=0; i<n_threads; i++)
        {
            unsigned int
                lower_bound = i*slice_size,
                upper_bound = (i+1)*slice_size;
            if (lower_bound == 0)
                lower_bound = 1;
            if (upper_bound > s.size())
                upper_bound = s.size();
            threads.push_back(std::thread(walsh_spectrum_cols_count,
                                          std::ref(local_counts[i]),
                                          s,
                                          lower_bound,
                                          upper_bound));

        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            for (auto &entry : local_counts[i])
                count[entry.first] += entry.second ;
        }
    }
    dict result;
    for (auto &entry : count)
        result[entry.first] = count[entry.first] ;
    return result;
}
