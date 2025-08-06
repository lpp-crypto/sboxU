#include "../sboxU.hpp"


Integer cpp_Spectrum::maximum() const
{
    if (content.empty())
        return 0;
    else
    {
        Integer result;
        for(auto &row : content)
            if ((row.first > result) && (row.second > 0))
                result = row.first;
        return result;
    }
}


Integer cpp_Spectrum::size() const
{
    return content.size();
}


std::map<Integer, Integer>::iterator cpp_Spectrum::iterator() 
{
    return content.begin();
}

    
Integer cpp_Spectrum::operator[] (const Integer x) 
{
    return content[x];
}

    
void cpp_Spectrum::operator+= (const cpp_Spectrum sp) 
{
    for(auto &row : sp.content)
        incr_by_amount(row.first, row.second);
}


std::vector<Integer> cpp_Spectrum::keys() const
{
    std::vector<Integer> result;
    result.reserve(content.size());
    for(auto row : content)
        result.push_back(row.first);
    return result;
}


void cpp_Spectrum::incr(const Integer entry)
{
    content[entry] ++;
}

    
void cpp_Spectrum::incr_by_amount(const Integer entry, const Integer amount)
{
    content[entry] += amount;
}


void cpp_Spectrum::incr_by_counting(const std::vector<Integer> vector_to_count)
{
    for(auto &c : vector_to_count)
        incr(c);
}
