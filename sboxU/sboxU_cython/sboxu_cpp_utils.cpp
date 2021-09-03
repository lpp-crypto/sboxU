#include "sboxu_cpp_utils.hpp"

BinWord oplus_cpp(BinWord x, BinWord y)
{
    return (x ^ y);
}

unsigned int hamming_weight_cpp(BinWord x)
{
    return __builtin_popcount(x);
}


BinWord parity_cpp(BinWord x)
{
    return __builtin_parity(x);
}

BinWord scal_prod_cpp(BinWord x, BinWord y)
{
    return parity_cpp(x & y);
}

std::vector<BinWord> component_cpp(BinWord a, std::vector<BinWord> f)
{
    std::vector<BinWord> result(f.size(), 0);
    for (BinWord x=0; x<f.size(); x++)
        result[x] = scal_prod_cpp(a, f[x]);
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
        if ((count[v] == 1) or (v >= s.size()))
            return false;
        else
            count[v] = 1;
    }
    return true;
}


Sbox inverse_cpp(Sbox s)
{
    Sbox result(s.size(), 0);
    for (unsigned int x=0; x<s.size(); x++)
        result[s[x]] = x;
    return result;
}



// !SECTION! Rank of sets of vectors


Integer rank_of_vector_set_cpp(std::vector<BinWord> l)
{
    Integer result = 0;
    for (unsigned int i=0; i<l.size(); i++)
    {
        if (l[i] > 0)
        {
            result ++;
            for (unsigned int j=i+1; j<l.size(); j++)
            {
                BinWord y = l[i] ^ l[j];
                if (y < l[j])
                    l[j] = y;
            }
        }
    }
    return result;
}




// !SECTION! Simple verifications

void check_length_cpp(Sbox s)
{
    unsigned long int pow_2 = 1;
    while (pow_2 < s.size())
        pow_2 <<= 1;
    if (pow_2 != s.size())
        throw NotSboxSized(s.size());
}
