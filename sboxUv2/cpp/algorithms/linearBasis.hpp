#ifndef _ALGOS_LINEAR_BASIS_
#define _ALGOS_LINEAR_BASIS_

#include "../sboxU.hpp"
#include "../core/include.hpp"


class cpp_Linear_basis
{
private:
    std::map<Integer, BinWord> basis;
public:
    cpp_Linear_basis() : basis() {} ;

    cpp_Linear_basis(const std::vector<BinWord> & l) ;

    bool add_to_span(BinWord x);

    bool is_in_span(BinWord x) const ;

    std::vector<BinWord> get_basis() const;

    inline Integer rank() const
    {
        return basis.size();
    }
    
    std::vector<BinWord> span() const;

    cpp_Linear_basis image_by(const cpp_BinLinearMap & L) const;


    inline std::map<Integer, BinWord>::const_iterator begin() const
    {
        return basis.cbegin();
    }
    
    inline std::map<Integer, BinWord>::const_iterator end() const
    {
        return basis.cend();
    }

    bool operator==(const cpp_Linear_basis & other_basis) const;
    
    bool operator<(const cpp_Linear_basis & other_basis) const;
};


std::vector<BinWord> cpp_complete_basis(
    const std::vector<BinWord> basis,
    const unsigned int n
    );


#endif
