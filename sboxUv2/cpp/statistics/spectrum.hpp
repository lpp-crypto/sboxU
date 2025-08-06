#ifndef _SPECTRUM_
#define _SPECTRUM_

#include "../sboxU.hpp"

class cpp_Spectrum
{
public:

    std::map<Integer, Integer> content;
    
    cpp_Spectrum() : content() {} ;

    Integer maximum() const;

    Integer size() const;
    
    std::map<Integer, Integer>::iterator iterator() ;
    
    Integer operator[] (const Integer x) ;
    
    void operator+= (const cpp_Spectrum sp) ;

    std::vector<Integer> keys() const ;

    void incr(const Integer entry);
    
    void incr_by_amount(const Integer entry, const Integer amount);

    void incr_by_counting(const std::vector<Integer> vector_to_count);
};

#endif
