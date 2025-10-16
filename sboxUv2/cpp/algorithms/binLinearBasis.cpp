#include "./binLinearBasis.hpp"



cpp_BinLinearBasis::cpp_BinLinearBasis(const std::vector<BinWord> & l) :
    basis()
{
    for(auto &v : l)
        add_to_span(v);
}


bool cpp_BinLinearBasis::add_to_span(BinWord x)
/** The content of basis is corresponds to binary vectors indexed by
 * their MSBs. At every instant, these vectors satisfy the following
 * properties:
 * - they are linearly independent
 * - if a vector has its most significant bit at position i,
 *   then v_i=0 for all other vector v in basis.
 *
 * This function relies heavily on the following facts:
 *
 *   a < (a ^ b)  <=>  MSB(b) \not\in Supp(a)
 *
 * An std::map is iterated through in order of ascending keys.
 *
 * This allows us to operate as follows.
 * 
 * 1. We first go through the elements with an MSB smaller than or
 * equal to that of x. We extract these from x, meaning that if x has
 * a 1 where the basis vector b has its MSB, then we replace x by
 * x^b. If at any point we get that x==0, then we stop: it was already
 * in the span of the basis. If a vector b has the same MSB as x, then
 * we still replace x by x^b, but then we recompute its MSB.
 *
 * 2. We then need to extract the resulting x from all vectors with a
 * greater MSB.
 */
{
    if (x == 0)
        return false;
    else
    {
        Integer m = cpp_msb(x);    
        for(auto b : basis)
        {
            BinWord y = x ^ b.second;
            if (y == 0)
                return false;  // x was already in the span, we do nothing
            else if ((b.first <= m) and (y < x))
                x = y;   // if b is smaller than x, then we extract it from x
            else if ((b.first > m))
                break;              // there is nothing left to extract from x
        }
        // at this point, x has to be non-zero, not in the span of the
        // basis, and with a new MSB: we keep it.
        m = cpp_msb(x);
        basis[m] = x;
        // then, we extract it from all bigger vectors
        for(auto b: basis)
            if (b.first > m)  // equality is only possible for x, so we don't care
            {
                BinWord y = x ^ b.second;
                if (y < b.second)
                    basis[b.first] = y;
            }
        return true;
    }
}


bool cpp_BinLinearBasis::is_in_span(BinWord x) const
{
    Integer m = cpp_msb(x);    
    for(auto b : basis)
    {
        BinWord y = x ^ b.second;
        if (y == 0)
            return true; 
        else if ((b.first <= m) and (y < x))
            x = y;   // if b is smaller than x, then we extract it from x
        else if ((b.first > m))
            return false;
    }
    return (x == 0);
}


std::vector<BinWord> cpp_BinLinearBasis::get_basis() const
{
    std::vector<BinWord> result;
    result.reserve(basis.size());
    for (auto b : basis)
        result.push_back(b.second);
    return result;
}


std::vector<BinWord> cpp_BinLinearBasis::span() const
{
    unsigned int total_size = 1 << basis.size();
    std::vector<BinWord>
        result(total_size, 0),
        vects = get_basis();
    for(BinWord mask=1; mask<total_size; mask++) // no need to modify the entry in 0
        result[mask] = cpp_linear_combination(vects, mask);
    return result;
}

std::vector<BinWord> cpp_complete_basis(
    const std::vector<BinWord> basis,
    const unsigned int n
    )
{
    cpp_BinLinearBasis lb(basis);
    std::vector<BinWord> result(basis.cbegin(), basis.cend());
    result.reserve(n);
    for(BinWord x=1; result.size() < n; x*=2)
    {
        if (lb.add_to_span(x))
            result.push_back(x);
    }
    return result;
}


cpp_BinLinearBasis cpp_BinLinearBasis::image_by(const cpp_BinLinearMap & L) const
{
    cpp_BinLinearBasis result;
    for (auto &b : basis)
        result.add_to_span(L(b.second));
    return result;
}



bool cpp_BinLinearBasis::operator==(const cpp_BinLinearBasis & other_basis) const
{
    std::vector<BinWord>
        b1 = get_basis(),
        b2 = other_basis.get_basis();
    if (b1.size() != b2.size())
        return false;
    else
    {
        for(unsigned int i=0; i<b1.size(); i++)
            if (b1[i] != b2[i])
                return false;
        return true;
    }
    // if (basis.size() != other_basis.rank())
    //     return false;
    // else
    // {
    //     for (auto & b : other_basis)
    //         if (! basis.contains(b.first))
    //             return false;
    //         else if (basis.at(b.first) != b.second)
    //             return false;
    // }
    // return true;
}


bool cpp_BinLinearBasis::operator<(const cpp_BinLinearBasis & other_basis) const
{

    std::vector<BinWord>
        b1 = get_basis(),
        b2 = other_basis.get_basis();
    for(unsigned int i=0; ((i < b1.size()) && (i < b2.size())); i++)
    {
        if (b1[i] > b2[i])
            return false;
        else if (b1[i] < b2[i])
            return true;
    }
    return (b1.size() < b2.size());
    // auto b1 = basis.begin();
    // auto b2 = other_basis.begin();
    // // the following works because std::map iterators iterate in
    // // ascending order of the key, i.e. of the MSBs here.
    // while ((b1 != basis.end()) && (b2 != other_basis.end()))
    // {
    //     // comapring values
    //     if (b1->second < b2->second)
    //         return true;
    //     else if (b1->second > b2->second)
    //         return false;
    //     // values are identical, we keep going
    //     else
    //     {
    //         b1 ++;
    //         b2 ++;
    //     }
    // }
    // // at this point, either b1 or b2 has reached its end
    // if (b2 != basis.end())     // if b2 is not finished
    //     return true;
    // else // either b1 is the only one finished (in which case b1 is
    //      // shorter, thus smaller), or they are equal, in which case
    //      // the strict inequality is also false
    //     return false;
}
