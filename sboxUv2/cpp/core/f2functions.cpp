#include "f2functions.hpp"




// !SECTION! Linear combinations and their ranks

std::vector<BinWord> cpp_transpose (const std::vector<BinWord> & l)
{
    BinWord biggest_msb = 0;
    for (auto &v: l)
    {
        BinWord m = cpp_msb(v);
        if (m > biggest_msb)
            biggest_msb = m;
    }
    biggest_msb ++;             // position start at 0, not 1
    std::vector<BinWord> result(biggest_msb, 0);
    for (BinWord i=0; i<biggest_msb; i++)
    {
        for (BinWord j=0; j<l.size(); j++)
            if ((l[j] >> i) & 1)
                result[i] |= 1 << j;
    }
    return result;
}


BinWord cpp_linear_combination (const std::vector<BinWord> & l, const BinWord mask)
{
    BinWord
        result = 0,
        bit_selector = 1;
    for(auto &x : l)
    {
        if (mask & bit_selector)
            result ^= x;
        bit_selector <<= 1;
    }
    return result;
}


Integer cpp_rank_of_vector_set(std::vector<BinWord> l)
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

std::vector<int> cpp_to_bin(const BinWord x, int n) {
    std::vector<int> bits;
    bits.reserve(n);
    for (int i =0; i <=n-1; i++) {
        bits.push_back((x >> i) & 1);
    }
    return bits;
}

BinWord cpp_from_bin(const std::vector<int>& v) {
    BinWord y = 0;
    for (int i =0; i < v.size(); i++) {
        y |= (static_cast<BinWord>(v[i] & 1) << i);
    }
    return y;
}

BinWord cpp_circ_shift(const BinWord x, int n, int shift){
    shift = ((shift % n) + n) % n;
    BinWord mask =  (1<<n) - 1;
    return (x >> shift) | ((x << (n - shift))& mask ); 
}

