#include "./binLinearBasis.hpp"
#include "./bigvectors.hpp"
#include "./linearSystem.hpp"


// !SECTION! cpp_BinLinearBasis


// !SUBSECTION! Simple methods
// ! The functions in this subsection correspond essentially to basic
// ! setters and getters; nothing deep is happening.

cpp_BinLinearBasis::cpp_BinLinearBasis(const std::vector<BinWord> & l) :
    basis()
{
    for(auto &v : l)
        add_to_span(v);
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

cpp_BinLinearBasis cpp_BinLinearBasis::image_by(const cpp_F2AffineMap & L) const
{
    if (L.is_linear())
    {
        cpp_BinLinearBasis result;
        for (auto &b : basis)
            result.add_to_span(L(b.second));
        return result;
    }
    else
    {
        throw std::runtime_error("A cpp_BinLinearBasis cannot be the input of a non-linear cpp_F2AffineMap");
    }
}


// !SUBSECTION! Combining bases

cpp_BinLinearBasis cpp_BinLinearBasis::operator+(const cpp_BinLinearBasis & L) const
{
    cpp_BinLinearBasis result;
    for(auto & b: basis)
        result.add_to_span(b.second);
    for(auto & b: L.get_basis())
        result.add_to_span(b);
    return result;
}



// !SUBSECTION! add_to_span: where the main difficulty lies
// ! add_to_span is the main method that modifies the state of a
// ! cpp_BinLinearBasis, it implements an algorithm that is
// ! non-trivial.

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

cpp_BinLinearBasis cpp_BinLinearBasis::intersection(const cpp_BinLinearBasis blb_1) const
{   std::vector <BinWord> basis_1=blb_1.get_basis();
    BinWord n1=blb_1.rank();
    BinWord n2=basis.size();
    cpp_F2LinearSystem system=cpp_F2LinearSystem(n1+n2,true);
    BinWord msb_1=blb_1.msb();
    BinWord msb_2=basis.rbegin()->first;
    for (BinWord i =0; i <= std::max(msb_1,msb_2); i++){
        std::vector<BinWord> tmp;
        BinWord ct=0;
        for (auto b : basis){
            if ((b.second>>i)&1){
                tmp.push_back(ct);
            }
            ct++;
        }
        for (auto b : basis_1){
            if ((b>>i)&1){
                tmp.push_back(ct);
            }
            ct++;
        }
        system.add_equation(tmp);
    }   
    std::vector<cpp_BigF2Vector> kernel=system.kernel();
    cpp_BinLinearBasis res=cpp_BinLinearBasis();
    for (auto v : kernel){// For each vector in the kernel, we add the corresponding vector in the intersection
        BinWord x = 0;
        BinWord ct=0;
        for (auto b : basis){
            if (v.is_set(ct)){
                x=cpp_oplus(x,b.second);
            }
            ct ++;
        }
        res.add_to_span(x);
    }
    return res;
}


// !SECTION! Helper functions 

std::vector<BinWord> cpp_complete_basis(
    const cpp_BinLinearBasis & basis,
    const unsigned int n
    )
{
    cpp_BinLinearBasis lb = basis;
    std::vector<BinWord> result = basis.get_basis();
    result.reserve(n);
    for(BinWord x=1; result.size() < n; x<<=1)
    {
        if (lb.add_to_span(x))
            result.push_back(x);
    }
    return result;
}


bool cpp_is_sum_full_rank(
    const cpp_BinLinearBasis & b1,
    const cpp_BinLinearBasis & b2
    )
{
    cpp_BinLinearBasis total = b1;
    for (auto & v : b2.get_basis())
        if (not total.add_to_span(v))
            return false;
    return true;
}
