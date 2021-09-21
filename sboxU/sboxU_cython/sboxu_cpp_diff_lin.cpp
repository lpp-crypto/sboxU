/* Time-stamp: <2021-09-21 11:07:53 lperrin>
 *
 * LICENSE
 */ 

#include "sboxu_cpp_diff_lin.hpp"


// !SECTION! Differential properties

// !SUBSECTION! Internal routines 


std::vector<Integer> ddt_row(const Sbox s, const BinWord a)
{
    std::vector<Integer> result(s.size(), 0);
    for (unsigned int x=0; x< s.size(); x++)
        result[s[x^a] ^ s[x]] ++ ;
    return result;
};


void ddt_rows_count(
    std::map<Integer,Integer> &result,
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

bool is_ddt_row_max_smaller_than(const Sbox s, const BinWord a, const Integer u)
{
    std::vector<Integer> row(s.size(), 0);
    for (unsigned int x=0; x<s.size(); x++)
    {
        BinWord d_out = s[x^a] ^ s[x];
        row[d_out] ++ ;
        if (row[d_out] > u)
            return false;
    }
    return true;
}


// !SUBSECTION! Python-facing functions 

std::vector< std::vector<Integer> > ddt_cpp(const Sbox s)
{
	std::vector< std::vector<Integer> > table_ddt ;
	check_length_cpp(s);
	for (unsigned int i = 0 ; i < s.size(); i++)
	{
		table_ddt.push_back(ddt_row(s,i));
	}
	return table_ddt ;
}

std::map<Integer,Integer> differential_spectrum_fast(const Sbox  s, const unsigned int n_threads)
{
    check_length_cpp(s);
    std::map<Integer,Integer> count;
    if (n_threads == 1)
    {
        // small S-Box
        ddt_rows_count(std::ref(count), s, 1, s.size());
    }
    else
    {
        std::vector<std::thread> threads;
        std::vector<std::map<Integer,Integer> > local_counts(n_threads);
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
    return count;
}


bool is_differential_uniformity_smaller_than_cpp(const Sbox s, const Integer u)
{
    check_length_cpp(s);
    for (unsigned int a=1; a<s.size(); a++)
        if (is_ddt_row_max_smaller_than(s, a, u) == false)
            return false;
    return true;
}


std::vector<std::map<Integer, Integer> > c_differential_spectra_cpp(
    const Sbox s,
    const Sbox l_table,
    const Sbox e_table)
{
    std::vector<std::map<Integer, Integer> > result(1, std::map<Integer, Integer>());
    result.push_back(differential_spectrum_fast(s, 1));
    unsigned int modulus = s.size() - 1;
    for (unsigned int c=2; c<s.size(); c++)
    {
        std::map<Integer, Integer> spectrum;
        Integer log_c = l_table[c];
        for(unsigned int a=0; a<s.size(); a++)
        {
            std::vector<Integer> c_ddt_row(s.size(), 0);
            for(unsigned int x=0; x<s.size(); x++)
            {
                Integer b = s[x ^ a] ^ e_table[(log_c + l_table[s[x]]) % modulus];
                c_ddt_row[b] += 1;
            }
            for(auto & val : c_ddt_row)
                spectrum[val] += 1;
        }
        result.push_back(spectrum);
    }
    return result;
}


// !SECTION! Linear properties

Sbox invert_lat_cpp(const std::vector< std::vector<Integer> > l, const unsigned int n)
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
                if (scal_prod_cpp(a, x) == 0)
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
            f[x] = scal_prod_cpp(b, s[x]);
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
        Sbox f(component_cpp(b, s)) ;
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
        Sbox f(component_cpp(b, s)) ;
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

std::vector< std::vector<Integer> > lat_cpp(const Sbox s)
{
    check_length_cpp(s);
    std::vector<std::vector<Integer> > table(s.size(), std::vector<Integer>(s.size(), 0));
    // generating table
    for (unsigned int b=0; b<s.size(); b++)
    {
        // computing one coordinate
        Sbox f(s.size(), 0);
        for (unsigned int x=0; x<f.size(); x++)
            f[x] = scal_prod_cpp(b, s[x]);
        // Walsh transform
        std::vector<Integer> w(walsh_spectrum_coord(f));
        table[b].assign(w.begin(), w.end()); 
    }
    std::vector<std::vector<Integer>> result ; 
    for (unsigned int a=0; a<s.size(); a++)
    {
        std::vector<Integer> row(s.size(), 0);
        for (unsigned int b=0; b<s.size(); b++)
            row[b] = table[b][a];
        result.push_back(row);
    }
    return result;
}

std::map<Integer, Integer> walsh_spectrum_fast_cpp(const Sbox s, const unsigned int n_threads)
{
    check_length_cpp(s);
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
    return count;
}

std::vector<BinWord> lat_zeroes_cpp(
    const Sbox s,
    const unsigned int n,
    const unsigned int n_threads)
{
    check_length_cpp(s);
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
    check_length_cpp(s);
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


// !SECTION! Quadratic functions


Sbox ortho_derivative_fast(const Sbox& s)
{
    check_length_cpp(s);
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
                if (scal_prod_cpp(ortho, b) == 0)
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

std::vector< std::vector<Integer>> bct_cpp(const Sbox s)
{
	std::vector< std::vector<Integer>> table_bct ;
	check_length_cpp(s) ;
	Sbox s_inv(inverse_cpp(s)) ;
	for(unsigned int a = 0; a < s.size(); a++)
		table_bct.push_back(bct_row(s,s_inv,a)) ;
	return table_bct ;
}

std::map<Integer, Integer> bct_spectrum_fast_cpp(const Sbox s, const unsigned int n_threads)
{
    check_length_cpp(s);
    Sbox s_inv(inverse_cpp(s)) ;
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
    return count;
}
