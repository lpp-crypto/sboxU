/* Time-stamp: <2018-04-18 17:51:01 lperrin>
 *
 * LICENSE
 */ 


#include "sboxu_cpp_utils.hpp"
using namespace boost::python;


uint32_t oplus_cpp(uint32_t x, uint32_t y)
{
    return (x ^ y);
}

unsigned int hamming_weight(uint32_t x)
{
    unsigned int result = 0;
    while (x != 0)
    {
        if (x & 1)
            result += 1;
        x >>= 1;
    }
    return result;
}


Integer scal_prod(uint32_t x, uint32_t y)
{
    Integer result = 0;
    uint32_t conjunction = x & y;
    while (conjunction != 0)
    {
        if (conjunction & 1)
            result ^= 1;
        conjunction >>= 1;
    }
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

list random_permutation(unsigned int n)
{
    return vec_2_lst_int(random_permutation_cpp(n)); 
}


bool is_permutation_cpp(Sbox s)
{
    std::map<uint32_t, Integer> count;
    for (auto &v : s)
    {
        if (count[v] == 1)
            return false;
        else
            count[v] = 1;
    }
    return true;
}


bool is_permutation(list s)
{
    return is_permutation_cpp(lst_2_vec_int(s));
}


Sbox inverse(Sbox s)
{
    Sbox result(s.size(), 0);
    for (unsigned int x=0; x<s.size(); x++)
        result[s[x]] = x;
    return result;
}


// !SECTION! Conversion C++/Python 

std::vector<Integer> lst_2_vec_int(const list &l)
{
    std::vector<Integer> result;
    for (unsigned int i=0; i<len(l); i++)
        result.push_back(extract<Integer>(l[i]));
    return result;
}


list vec_2_lst_int(const std::vector<Integer> l)
{
    list result;
    for (auto &v : l)
        result.append<Integer>(v);
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

