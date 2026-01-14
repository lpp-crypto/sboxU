#ifndef _ALGOS_LINEAR_BASIS_
#define _ALGOS_LINEAR_BASIS_

#include "../common.hpp"
#include "../core/include.hpp"


class cpp_BinLinearBasis
{
private:
    std::map<Integer, BinWord> basis;
public:
    cpp_BinLinearBasis() : basis() {} ;

    cpp_BinLinearBasis(const std::vector<BinWord> & l) ;

    // !TODO! destructor for BinLinearBasis 
    // inline void destruct()
    // {
    //     for (auto it : basis)
    //         delete basis.second;
    // }

    bool add_to_span(BinWord x);

    bool is_in_span(BinWord x) const ;

    std::vector<BinWord> get_basis() const;

    inline Integer rank() const
    {
        return basis.size();
    }
    
    std::vector<BinWord> span() const;

    cpp_BinLinearBasis image_by(const cpp_BinLinearMap & L) const;


    inline std::map<Integer, BinWord>::const_iterator begin() const
    {
        return basis.cbegin();
    }
    
    inline std::map<Integer, BinWord>::const_iterator end() const
    {
        return basis.cend();
    }


    // !SUBSECTION! Boolean operators 

    inline bool operator==(const cpp_BinLinearBasis & other_basis) const
    {
        if (basis.size() != other_basis.rank())
            return false;
        else
        {
            for (auto & b : other_basis)
                if (! basis.contains(b.first))
                    return false;
                else if (basis.at(b.first) != b.second)
                    return false;
        }
        return true;
    }


    inline bool operator<(const cpp_BinLinearBasis & other_basis) const
    {
        auto b1 = basis.begin();
        auto b2 = other_basis.begin();
        // the following works because std::map iterators iterate in
        // ascending order of the key, i.e. of the MSBs here.
        while ((b1 != basis.end()) && (b2 != other_basis.end()))
        {
            // values are identical, we keep going
            if (b1->second == b2->second)
            {
                b1 ++;
                b2 ++;
            }
            else
                // comparing values
                return (b1->second < b2->second);
        }
        // at this point, either b1 or b2 has reached its end. If b1 is
        // the only one finished, then b1 is shorter, thus smaller. If b2
        // has finished, then one of two things is happening.  b2 has
        // finished, meaning b1 is longer, and thus larger, so we return
        // false. Alternatively, they are of equal length, and thus at
        // this point equal: we also return false (we test for strict
        // inequality.
        return (b2 != other_basis.end());
    }
};


std::vector<BinWord> cpp_complete_basis(
    const cpp_BinLinearBasis & basis,
    const unsigned int n
    );


#endif
