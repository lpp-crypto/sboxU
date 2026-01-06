#ifndef _BIG_VECTOR_
#define _BIG_VECTOR_

#include "../sboxU.hpp"

class cpp_SparseVector
{
public:
    std::vector<unsigned int> indices;
    
    cpp_SparseVector(std::vector<unsigned int> _indices) :
        indices(_indices.begin(), _indices.end())
    {
        if (indices.size() > 0)
            std::sort(indices.begin(), indices.end(), std::greater{});
    };


    bool is_zero() const
    {
        return (indices.size() > 0);
    };
    
    BinWord operator[](const unsigned int index) const
    {
        if (std::binary_search(indices.begin(), indices.end(), index))
            return 1;
        else
            return 0;
    }

    
    BinWord operator^(const cpp_SparseVector & other_vector) const
    {
        std::vector<unsigned int> result_indices;
        result_indices.reserve(hamming_weight() + other_vector.hamming_weight());
        for(auto & index : indices)
            if (other_vector[index] == 0)
                result_indices.push_back(index);
        for(auto & index : other_vector.indices)
            if (indices[index] == 0)
                result_indices.push_back(index);
        return cpp_SparseVector(result_indices);
    }
};


bool operator<(const cpp_SparseVector & v1,
               const cpp_SparseVector & v2) const
{
    return std::lexicographical_compare(
        v1.indices.begin(), v1.indices.end(),
        v2.indices.begin(), v2.indices.end()
        );
}


bool operator==(const cpp_SparseVector & v1,
                const cpp_SparseVector & v2) const
{
    return (v1.indices == v2.indices);
}


#endif
