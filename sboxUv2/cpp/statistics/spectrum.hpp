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

    inline Integer size() const
    {
        return content.size();
    }
    
    inline std::map<Integer, Integer>::const_iterator cbegin() const
    {
        return content.cbegin();
    };

    inline std::map<Integer, Integer>::const_iterator cend() const
    {
        return content.cend();
    };
    
    inline Integer operator[] (const Integer x)
    {
        return content[x];
    };
    
    
    void operator+= (const cpp_Spectrum sp) ;

    std::vector<Integer> keys() const ;

    inline void incr(const Integer entry)
    {
        content[entry] ++;
    }
    
    inline void incr_by_amount(const Integer entry, const Integer amount)
    {
        content[entry] += amount;
    }

    void incr_by_counting(const std::vector<Integer> & vector_to_count);
};

#endif
