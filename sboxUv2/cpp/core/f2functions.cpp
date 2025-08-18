#include "f2functions.hpp"




// !SECTION! Linear combinations and their ranks

BinWord cpp_linear_combination (const std::vector<BinWord> & l, BinWord mask)
{
    BinWord result = 0;
    for(auto &x : l)
    {
        if (mask & 1)
            result ^= x;
        mask >>= 1;
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


