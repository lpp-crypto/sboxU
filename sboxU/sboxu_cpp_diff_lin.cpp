/* Time-stamp: <2018-05-16 14:45:06 lperrin>
 *
 * LICENSE
 */ 

#include "sboxu_cpp_diff_lin.hpp"
using namespace boost::python;


// !SECTION! Differential properties

// !SUBSECTION! Internal routines 

std::vector<Integer> ddt_row(const Sbox s, const BinWord a)
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
    const BinWord a_min,
    const BinWord a_max)
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
    Sbox s(lst_2_vec_BinWord(l)) ;
    check_length(s);
    for (unsigned int a=0; a<s.size(); a++)
        result.append(vec_2_lst_Integer(ddt_row(s, a)));

    return result;
}


dict differential_spectrum_fast(const list& l, const unsigned int n_threads)
{
    Sbox s(lst_2_vec_BinWord(l)) ;
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


Sbox invert_lat_cpp(const std::vector<std::vector<Integer> > l, const unsigned int n)
{
    Sbox result(l.size(), 0);
    for (unsigned int i=0; i<n; i++)
    {
        Integer b = (1 << i), sum;
        for (unsigned int x=0; x<l.size(); x++)
        {
            sum = 0;
            for (unsigned int a=0; a<l.size(); a++)
            {
                if (scal_prod(a, x) == 0)
                    sum += l[a][b];
                else
                    sum -= l[a][b];
            }
            if (sum < 0)
                result[x] |= (1 << i) ;
        }
    }
    return result;
}


// !SUBSECTION! Internal routines 

/* The algorithm (and its notations) come from:
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
    const BinWord b_min,
    const BinWord b_max)
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

void lat_zeroes_in_columns(
    std::vector<BinWord> &result,
    const Sbox s,
    const BinWord b_min,
    const BinWord b_max,
    const unsigned int n)
{
    for (unsigned int b=b_min; b<b_max; b++)
    {
        // computing one coordinate
        Sbox f(component(b, s)) ;
        // Walsh transform
        std::vector<Integer> w(walsh_spectrum_coord(f));
        // getting zeroes
        for (unsigned int a=0; a<w.size(); a++)
            if (w[a] == 0)
                result.push_back((a << n) | b) ;
    }
}




// !SUBSECTION! Python-facing functions 

list lat(const list& l)
{
    Sbox s(lst_2_vec_BinWord(l)) ;
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
        result.append<list>(vec_2_lst_Integer(row));
    }
    return result;
}


dict walsh_spectrum_fast(const list& l, const unsigned int n_threads)
{
    Sbox s(lst_2_vec_BinWord(l)) ;
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


list lat_zeroes_fast(const list& l,
                     const unsigned int n,
                     const unsigned int n_threads)
{
    Sbox s(lst_2_vec_BinWord(l)) ;
    check_length(s);
    std::vector<BinWord> zeroes;
    if (n_threads == 1)
    {
        // small S-Box
        lat_zeroes_in_columns(std::ref(zeroes), s, 0, s.size(), n);
    }
    else
    {
        std::vector<std::thread> threads;
        std::vector<std::vector<BinWord> > local_zeroes(n_threads);
        unsigned int slice_size = s.size()/n_threads;
        for (unsigned int i=0; i<n_threads; i++)
        {
            unsigned int
                lower_bound = i*slice_size,
                upper_bound = (i+1)*slice_size;
            if (upper_bound > s.size())
                upper_bound = s.size();
            threads.push_back(std::thread(lat_zeroes_in_columns,
                                          std::ref(local_zeroes[i]),
                                          s,
                                          lower_bound,
                                          upper_bound,
                                          n));

        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            zeroes.insert(zeroes.end(),
                          local_zeroes[i].begin(),
                          local_zeroes[i].end());
        }
    }
    return vec_2_lst_BinWord(zeroes) ;
}
