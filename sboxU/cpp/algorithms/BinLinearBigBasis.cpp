#include "./BinLinearBigBasis.hpp"
#include "./bigvectors.hpp"


cpp_BinLinearBigBasis::cpp_BinLinearBigBasis(const std::vector<cpp_BigF2Vector> & l) :
    basis()
{
    if (l.empty())
    {
        dimension = 0;
        return;
    }
    dimension = l[0].size();
    for(const auto &v : l){
        if (v.size() != dimension)
            throw std::invalid_argument(
                "cpp_BinLinearBigBasis: inconsistent vector dimensions"
            );
        add_to_span(v);
}
    }
        

cpp_BinLinearBigBasis::cpp_BinLinearBigBasis(const std::vector<std::vector<BoolBlock>> & l, unsigned int n) :
    basis(),dimension(n)
{
    for (const auto &v : l)
    {
        if (v.size() * BLOCK_SIZE < dimension)
            throw std::invalid_argument(
                "cpp_BinLinearBigBasis: insufficient blocks for dimension"
            );

        add_to_span(v);
    }
}


bool cpp_BinLinearBigBasis::add_to_span(cpp_BigF2Vector big_x)
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
{   if (big_x.size() != dimension)
        throw std::invalid_argument(
            "add_to_span: vector dimension mismatch"
        );

    if (big_x.is_zero())
        return false;
    else
    {
        Integer m = big_x.get_msb();    
        for(auto b : basis)
        {
            cpp_BigF2Vector big_y = big_x ^ b.second;
            if (big_y.is_zero())
                return false;  // x was already in the span, we do nothing
            else if ((b.first <= m) and (big_y < big_x))
                big_x = big_y;   // if b is smaller than x, then we extract it from x
            else if ((b.first > m))
                break;              // there is nothing left to extract from x
        }
        // at this point, x has to be non-zero, not in the span of the
        // basis, and with a new MSB: we keep it.
        m = big_x.get_msb();
        basis[m] = big_x;
        // then, we extract it from all bigger vectors
        for(auto b: basis)
            if (b.first > m)  // equality is only possible for x, so we don't care
            {
                cpp_BigF2Vector big_y = big_x ^ b.second;
                if (big_y < b.second)
                    basis[b.first] = big_y;
            }
        return true;
    }
}



bool cpp_BinLinearBigBasis::add_to_span(std::vector<BoolBlock> x)
{   
    cpp_BigF2Vector big_x(dimension);

    for (unsigned int i = 0; i < x.size(); ++i)
        if (x[i])
            big_x.set_to_1(i);
    return add_to_span(big_x);
}


bool cpp_BinLinearBigBasis::is_in_span(std::vector<BoolBlock> x) const
{   
    cpp_BigF2Vector big_x=cpp_BigF2Vector(x,x.size()*BLOCK_SIZE);
    Integer m = big_x.get_msb();    
    for(auto b : basis)
    {
        cpp_BigF2Vector big_y = big_x ^ b.second;
        if (big_y.is_zero())
            return true; 
        else if ((b.first <= m) and (big_y < big_x))
            big_x = big_y;   // if b is smaller than x, then we extract it from x
        else if ((b.first > m))
            return false;
    }
    return (big_x.is_zero());
}


std::vector<Bytearray> cpp_BinLinearBigBasis::get_basis() const
{
    std::vector<Bytearray> result;
    result.reserve(basis.size());
    for (auto b : basis)
        result.push_back((b.second).to_bits_merlin());
    return result;
}


// std::vector<cpp_BigF2Vector> cpp_BinLinearBigBasis::span() const
// {
//     unsigned int total_size = 1 << basis.size();
//     std::vector<cpp_BigF2Vector>
//         result(total_size, 0),
//         vects = get_basis();
//     for(cpp_BigF2Vector mask=1; mask<total_size; mask++) // no need to modify the entry in 0
//         result[mask] = cpp_linear_combination(vects, mask);
//     return result;
// }

// std::vector<cpp_BigF2Vector> cpp_complete_basis(
//     const cpp_BinLinearBigBasis & basis,
//     const unsigned int n
//     )
// {
//     cpp_BinLinearBigBasis lb = basis;
//     std::vector<cpp_BigF2Vector> result = basis.get_basis();
//     result.reserve(n);
//     for(cpp_BigF2Vector x=1; result.size() < n; x*=2)
//     {
//         if (lb.add_to_span(x))
//             result.push_back(x);
//     }
//     return result;
// }


// cpp_BinLinearBigBasis cpp_BinLinearBigBasis::image_by(const cpp_BinLinearMap & L) const
// {
//     cpp_BinLinearBigBasis result;
//     for (auto &b : basis)
//         result.add_to_span(L(b.second));
//     return result;
// }

