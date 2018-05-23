/* Time-stamp: <2018-05-16 14:40:43 lperrin>
 *
 * LICENSE
 */ 


#include "sboxu_cpp_utils.hpp"
using namespace boost::python;


BinWord oplus_cpp(BinWord x, BinWord y)
{
    return (x ^ y);
}

unsigned int hamming_weight(BinWord x)
{
    return __builtin_popcount(x);
}


BinWord parity(BinWord x)
{
    return __builtin_parity(x);
}

BinWord scal_prod(BinWord x, BinWord y)
{
    return parity(x & y);
}

std::vector<BinWord> component(BinWord a, std::vector<BinWord> f)
{
    std::vector<BinWord> result(f.size(), 0);
    for (BinWord x=0; x<f.size(); x++)
        result[x] = scal_prod(a, f[x]);
    return result;
}


// !SECTION! Generating and studying permutations 

Sbox random_permutation_cpp(unsigned int n)
{
    // Initializing PRNG
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    // Initializing LUT 
    Integer size = 1;
    for (unsigned int i=0; i<n; i++)
        size *= 2;
    Sbox result(size, 0);
    // Knuth shuffle
    unsigned int j = 0;
    for (unsigned int i=0; i<result.size(); i++)
    {
        result[i] = i;
        j = generator() % (i+1);
        std::swap(result[i], result[j]);
    }
    return result;        
}


bool is_permutation_cpp(Sbox s)
{
    std::map<BinWord, Integer> count;
    for (auto &v : s)
    {
        if (count[v] == 1)
            return false;
        else
            count[v] = 1;
    }
    return true;
}


Sbox inverse(Sbox s)
{
    Sbox result(s.size(), 0);
    for (unsigned int x=0; x<s.size(); x++)
        result[s[x]] = x;
    return result;
}


// !SECTION! Conversion C++/Python 

std::vector<Integer> lst_2_vec_Integer(const list &l)
{
    std::vector<Integer> result;
    for (unsigned int i=0; i<len(l); i++)
        result.push_back(extract<Integer>(l[i]));
    return result;
}


std::vector<BinWord> lst_2_vec_BinWord(const list &l)
{
    std::vector<BinWord> result;
    for (unsigned int i=0; i<len(l); i++)
        result.push_back(extract<BinWord>(l[i]));
    return result;
}


list vec_2_lst_Integer(const std::vector<Integer> l)
{
    list result;
    for (auto &v : l)
        result.append<Integer>(v);
    return result;
}

list vec_2_lst_BinWord(const std::vector<BinWord> l)
{
    list result;
    for (auto &v : l)
        result.append<BinWord>(v);
    return result;
}


// !SECTION! Simple verifications

void check_length(Sbox s)
{
    unsigned long int pow_2 = 1;
    while (pow_2 < s.size())
        pow_2 <<= 1;
    if (pow_2 != s.size())
        throw NotSboxSized(s.size());
}

