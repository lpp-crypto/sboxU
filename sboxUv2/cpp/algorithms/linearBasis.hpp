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

    void add_to_span(BinWord x);

    bool is_in_span(BinWord x) const ;

    std::vector<BinWord> get_basis() const;

    inline Integer rank() const
    {
        return basis.size();
    }
    
    std::vector<BinWord> span() const;
};


#endif
