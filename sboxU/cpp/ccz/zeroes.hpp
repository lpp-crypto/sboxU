#ifndef _CCZ_ZEROES_
#define _CCZ_ZEROES_

#include "../common.hpp"
#include "../core/include.hpp"
#include "../statistics/include.hpp"
#include "../algorithms/include.hpp"



class cpp_WalshZeroesSpaces
{
private:
    BinWord mask;
    BinWord n;
    unsigned int total_size;
public:
    std::vector<cpp_BinLinearBasis> bases;
    std::vector<cpp_F2AffineMap> mappings;

    cpp_WalshZeroesSpaces() : bases(0), mappings(0) {} ;
    
    cpp_WalshZeroesSpaces(
        const std::vector<cpp_BinLinearBasis> & _bases,
        const unsigned int _n,
        const unsigned int _total_size
        ):
        mask((1 << _n) - 1),
        n(_n),
        total_size(_total_size),
        bases(_bases.begin(), _bases.end()),
        mappings(0)
    {};
    
    cpp_WalshZeroesSpaces(
        const cpp_S_box & s,
        const unsigned int n_threads
        );

    void destruct()
    {
        bases.clear();
        bases.shrink_to_fit();
        mappings.clear();
        mappings.shrink_to_fit();
    }

    ~cpp_WalshZeroesSpaces()
    {
        destruct();
    }

    void init_mappings();

    
    /** @brief Computes the set of the admissible mappings corresponding to these Walsh zeroes and reduces it using the set of automorphisms provided.
     *
     * At the end of its execution, the this->mappings attribute is properly initialised. If the set of automorphisms provided is complete, then there is a bijection between the mappings contained in the this->mappings attribute and the EA-classes contained in the CCZ-class of the function assoicated to this instance.
     *
     * @param automorphisms An std::vector of affine mappings that is expected to contain all the automorphisms of the function associated to this cpp_WalshZeroesSpaces instance. This is not verified!
     */
    void init_mappings(const std::vector<cpp_F2AffineMap> & automorphisms);

    
    /** @brief Computes the set of the admissible mappings corresponding to these Walsh zeroes and reduces it using the sets of automorphisms provided (assuming that the complete set of the automorphisms is the product of the two sets given).
     *
     * At the end of its execution, the this->mappings attribute is properly initialised. If the set of automorphisms provided is complete, then there is a bijection between the mappings contained in the this->mappings attribute and the EA-classes contained in the CCZ-class of the function assoicated to this instance.
     *
     * @param automorphisms_1 An std::vector of affine mappings that is expected to contain automorphisms of the function associated to this cpp_WalshZeroesSpaces instance. This is not verified!
     *
     * @param automorphisms_2 Another std::vector of affine mappings that is also expected to contain automorphisms of the function associated to this cpp_WalshZeroesSpaces instance. Any autmorphism of the function should be written as A_1*A_2, where "*" is the composition, A_1 is in automorphisms_1, and A_2 is in automorphisms_2. None of this is verified!
     */
    void init_mappings(
        const std::vector<cpp_F2AffineMap> & automorphisms_1,
        const std::vector<cpp_F2AffineMap> & automorphisms_2
        );

    cpp_WalshZeroesSpaces image_by(const cpp_F2AffineMap & L) const;
    
    cpp_Spectrum thickness_spectrum() const;
};


cpp_Spectrum cpp_thickness_spectrum(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

#endif
