#include "./linearBasis.hpp"



cpp_Linear_basis::cpp_Linear_basis(const std::vector<BinWord> & l) :
    basis()
{
    for(auto &v : l)
        add_to_span(v);
}


void cpp_Linear_basis::add_to_span(BinWord x)
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
    Integer m = cpp_msb(x);    
    for(auto b : basis)
    {
        BinWord y = x ^ b.second;
        if (y == 0)
            return;  // x was already in the span, we do nothing
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
}

bool cpp_Linear_basis::is_in_span(BinWord x) const
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


std::vector<BinWord> cpp_Linear_basis::get_basis() const
{
    std::vector<BinWord> result;
    result.reserve(basis.size());
    for (auto b : basis)
        result.push_back(b.second);
    return result;
}


std::vector<BinWord> cpp_Linear_basis::span() const
{
    unsigned int total_size = 1 << basis.size();
    std::vector<BinWord>
        result(total_size, 0),
        vects = get_basis();
    for(BinWord mask=1; mask<total_size; mask++) // no need to modify the entry in 0
        result[mask] = cpp_linear_combination(vects, mask);
    return result;
}
