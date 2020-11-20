/* Time-stamp: <2020-11-20 13:10:22 leo>
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
    result.reserve(result.size() + b_max - b_min) ;
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


void do_lat_columns_contain_zero(
    std::vector<bool> &result,
    const Sbox s,
    const BinWord b_min,
    const BinWord b_max)
{
    result.reserve(result.size() + b_max - b_min) ;
    for (unsigned int b=b_min; b<b_max; b++)
    {
        // computing one coordinate
        Sbox f(component(b, s)) ;
        // Walsh transform
        std::vector<Integer> w(walsh_spectrum_coord(f));
        bool contains_zero = false;
        for (auto & c : w)
            if (c == 0)
            {
                contains_zero = true;
                break;
            }
        result.push_back(contains_zero) ;
    }
}


// !SUBSECTION! High level parallelized functions 

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


std::vector<BinWord> lat_zeroes_cpp(
    const Sbox s,
    const unsigned int n,
    const unsigned int n_threads)
{
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
    return zeroes ;
}


std::vector<BinWord> projected_lat_zeroes_cpp(
    const Sbox s,
    const unsigned int n_threads)
{
    check_length(s);
    std::vector<bool> projection;
    if (n_threads == 1)
    {
        // small S-Box
        do_lat_columns_contain_zero(std::ref(projection), s, 0, s.size());
    }
    else
    {
        std::vector<std::thread> threads;
        std::vector<std::vector<bool> > local_projections(n_threads);
        unsigned int slice_size = s.size()/n_threads;
        for (unsigned int i=0; i<n_threads; i++)
        {
            unsigned int
                lower_bound = i*slice_size,
                upper_bound = (i+1)*slice_size;
            if (upper_bound > s.size())
                upper_bound = s.size();
            threads.push_back(std::thread(do_lat_columns_contain_zero,
                                          std::ref(local_projections[i]),
                                          s,
                                          lower_bound,
                                          upper_bound));

        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            projection.insert(projection.end(),
                          local_projections[i].begin(),
                          local_projections[i].end());
        }
    }
    std::vector<BinWord> projected_zeroes;
    for (unsigned int i=1; i<projection.size(); i++) // we purposefully leave out the zero
        if (projection[i])
            projected_zeroes.push_back(i) ;
    return projected_zeroes ;
}


// !SECTION! Boomerang properties

// !SUBSECTION! Internal routines 


// def bct(s):
//     n = int(log(len(s), 2))
//     s_inv = inverse(s)
//     table = [[0 for b in xrange(0, 2**n)] for a in xrange(0, 2**n)]
//     for b in xrange(0, 2**n):
//         T = [[] for y in xrange(0, 2**n)]
//         for x in xrange(0, 2**n):
//             y = oplus(x, s_inv[oplus(s[x], b)])
//             T[y].append(x)
//         for y in xrange(0, 2**n):
//             for x_i, x_j in itertools.product(T[y], T[y]):
//                 table[oplus(x_i, x_j)][b] += 1
//     return table


std::vector<Integer> bct_row(const Sbox s, const Sbox s_inv, const BinWord a)
{
    std::vector<Integer> result(s.size(), 0);
    std::vector<std::vector<BinWord> > xor_list(s.size(), std::vector<BinWord>(0));
    for (BinWord x=0; x<s.size(); x++)
    {
        BinWord y = x ^ s[s_inv[x] ^ a];
        xor_list[y].push_back(x);
    }
    for (BinWord y=0; y<s.size(); y++)
        for (BinWord &x1 : xor_list[y])
            for (BinWord &x2 : xor_list[y])
                result[x1 ^ x2] ++ ;
    return result;
}


void bct_rows_count(
    std::map<Integer, Integer> &result,
    const Sbox s,
    const Sbox s_inv,
    const BinWord a_min,
    const BinWord a_max)
{
    for (unsigned int a=a_min; a<a_max; a++)
    {
        std::vector<Integer> row(bct_row(s, s_inv, a));
        for (unsigned int i=1; i<row.size(); i++)
            result[row[i]] ++;
    }
}


// !SUBSECTION! Python-facing functions 

list bct(const list& l)
{
    list result;    
    Sbox s(lst_2_vec_BinWord(l)), s_inv(inverse(s)) ;
    check_length(s);
    for (unsigned int a=0; a<s.size(); a++)
        result.append(vec_2_lst_Integer(bct_row(s, s_inv, a)));

    return result;
}


dict bct_spectrum_fast(const list& l, const unsigned int n_threads)
{
    Sbox s(lst_2_vec_BinWord(l)), s_inv(inverse(s)) ;
    check_length(s);
    std::map<Integer, Integer> count;
    if (n_threads == 1)
    {
        // small S-Box
        bct_rows_count(std::ref(count), s, s_inv, 1, s.size());
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
            threads.push_back(std::thread(bct_rows_count,
                                          std::ref(local_counts[i]),
                                          s,
                                          s_inv,
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


// !SECTION! Quadratic functions


Sbox ortho_derivative_fast(const Sbox& s)
{
    Sbox result(s.size(), 0);
    for (unsigned int a=1; a<s.size(); a++)
    {
        // getting the hyperplane
        std::vector<Integer> row(ddt_row(s, a));
        std::vector<Integer> hyperplane;
        hyperplane.reserve(s.size() / 2);
        for(unsigned int b=1; b<s.size(); b++)
            if (row[b] != row[0])
                hyperplane.push_back(b);

        // we return an empty list if the function is not APN, which
        // is equivalent to all rows having exactly half of their
        // elements be non-zero
        if (hyperplane.size() < (s.size()/2))
            return Sbox(0);

        // bruteforcing "ortho" until it is orthogonal to all elements
        // in the hyperplane
        BinWord ortho = 1;
        bool found = false;
        while ((not found) and (ortho < s.size()))
        {
            found = true;
            for(auto &b : hyperplane)
                if (scal_prod(ortho, b) == 0)
                {
                    found = false;
                    break;
                }
            if (not found)
                ortho += 1;
        }
        // if we couldn't find an element orthogonal to the
        // hyperplane, then it is not a hyperplane and we return an
        // empty function
        if (found)
            result[a] = ortho;
        else
            return Sbox(0);
    }
    return result;
}


list ortho_derivative(const list& l)
{
    list result;    
    Sbox s(lst_2_vec_BinWord(l));
    check_length(s);
    return vec_2_lst_BinWord(ortho_derivative_fast(s));
}
