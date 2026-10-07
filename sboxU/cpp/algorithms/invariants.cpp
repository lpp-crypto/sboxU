#include "./invariants.hpp"
#include "bigvectors.hpp"
#include "BinLinearBigBasis.hpp"
#include "../core/s_box.hpp"
#include <bit>
#include <tuple>

std::vector<cpp_BigF2Vector> cpp_basis_invariants_from_cycles(const cpp_S_box S){
    
    std::vector<std::vector<BinWord>> cycles = cpp_cycle_decomposition(S);
    std::vector<cpp_BigF2Vector> result;
    int all_even = 1;
    unsigned int N = S.input_space_size();
    for (auto cycle : cycles){
        if ((cycle.size() % 2)==1){
            all_even=0;
            break;
        }
    }

    cpp_BigF2Vector big_vector;
    for (auto cycle : cycles) // The basis always contains invariants that are constant on the cycles of S.
    {
        big_vector =cpp_BigF2Vector(N);
        for (auto a : cycle){
            big_vector.set_to_1(a);
        }
        result.push_back(big_vector);
    }

    if (all_even)
    {
        big_vector = cpp_BigF2Vector(N);
        for (auto cycle : cycles)
            for (unsigned int i = 1; i < cycle.size(); i += 2)
                big_vector.set_to_1(cycle[i]);
        result.push_back(big_vector);
    }
    return result;
}

std::vector<BinWord> cpp_vectors_of_hamming_weight(BinWord h, BinWord n) // Generates the list of all vectors of hamming weight h and length n
{
    std::vector<BinWord> result; 
    if (h <0 || h> n) return result; 
    if (h==0) return {0};
    BinWord limit = 1 << n;
    BinWord x = (1<< h) -1;
    while (x < limit){
        result.push_back(x);
        BinWord c = x & -x;
        BinWord r = x + c;
        x = (((r ^ x) >> 2) / c) | r;
    }
    return result;

}

std::tuple<Lut,BinWord> yann_permutation(BinWord d, BinWord n) // Computes a permutation of [0,2^{n}-1] such that all intgers of hamming_weigth <= d appear first
{
    if (d >=n){
        throw std::runtime_error("In yann_permutation we need d<n");
    }
    else{
        const BinWord limit = 1 << n;
        Lut perm;
        perm.reserve(limit);
        for (BinWord h = 0; h <= d; ++h)
            for (auto x : cpp_vectors_of_hamming_weight(h, n))
                perm.push_back(x);
        BinWord bound = perm.size();
        for (BinWord x = 0; x < limit; ++x)
            if (std::popcount(x) > d)
                {perm.push_back(x);}
        return {std::move(perm),bound};
    }
}

std::vector<cpp_BigF2Vector> cpp_all_invariants_up_to_degree(const cpp_S_box S, BinWord d){
    std::vector < cpp_BigF2Vector> result;
    BinWord n = S.get_input_length();
    auto [perm,bound]= yann_permutation(d,n);
    Lut inv_perm= cpp_inverse(perm);
    cpp_BinLinearBigBasis basis = cpp_BinLinearBigBasis(1 << n);
    for (auto b : cpp_basis_invariants_from_cycles(S)){
        basis.add_to_span(apply_perm_BigF2Vector(mobius_transform(b,n),perm));
    }
    for (auto b : basis.get_basis())
    {
        if (b.get_msb()>=bound){
            break;
        }
        result.push_back(mobius_transform(apply_perm_BigF2Vector(b,inv_perm),n));
    }
    return result;
}