#include "spectrum.hpp"


Integer cpp_Spectrum::maximum() const
{
    if (content.empty())
        return 0;
    else
    {
        Integer result = 0;
        for(auto &row : content)
            if ((row.first > result) && (row.second > 0))
                result = row.first;
        return result;
    }
}
    
void cpp_Spectrum::operator+= (const cpp_Spectrum sp) 
{
    for(auto row=sp.cbegin(); row!=sp.cend(); row++)
        incr_by_amount(row->first, row->second);
}


std::vector<Integer> cpp_Spectrum::keys() const
{
    std::vector<Integer> result;
    result.reserve(content.size());
    for(auto row : content)
        result.push_back(row.first);
    return result;
}


void cpp_Spectrum::incr_by_counting(const std::vector<Integer> & vector_to_count)
{
    for(auto &c : vector_to_count)
        incr(c);
}


std::string cpp_Spectrum::content_string_repr() const
{
    std::stringstream result;
    result << "{";
    for(auto row : content)
        result << row.first << ":" << row.second << ",";
    result << "}";
    return result.str();
}


cpp_Spectrum cpp_Spectrum::absolute() const
{
    cpp_Spectrum result;
    for(auto row : content)
        result.incr_by_amount(llabs(row.first), row.second);
    return result;
}
