#ifndef _STAT_SPECTRUM_
#define _STAT_SPECTRUM_

#include "../sboxU.hpp"
#include "../s_box.hpp"


class cpp_Spectrum
{
private:
    std::map<Integer, Integer> content;
    
public:
    
    cpp_Spectrum() : content() {} ;

    Integer maximum() const;

    Integer size() const;
    
    inline std::map<Integer, Integer>::const_iterator cbegin() const
    {
        return content.cbegin();
    };

    inline std::map<Integer, Integer>::const_iterator cend() const
    {
        return content.cend();
    };
    
    Integer operator[] (const Integer x) ;
    
    void operator+= (const cpp_Spectrum sp) ;

    std::vector<Integer> keys() const ;

    void incr(const Integer entry);
    
    void incr_by_amount(const Integer entry, const Integer amount);

    void incr_by_counting(const std::vector<Integer> vector_to_count);
};

#endif
