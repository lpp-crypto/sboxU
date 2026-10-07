#include "./invariants.hpp"
#include "bigvectors.hpp"
#include "BinLinearBigBasis.hpp"
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

    if (all_even){ // If all cycles are of even length, the basis also contains invariants that are alternating on the cycles of S.
        for (auto cycle : cycles)
        {
            big_vector = cpp_BigF2Vector(N);

            for (int i=1; i < cycle.size(); i+=2)
            {
                big_vector.set_to_1(cycle[i]);
            }
            result.push_back(big_vector);
        }
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
    if (d > n){
        throw std::runtime_error("In yann_permutation we need d<n");
    }
    else{
        const BinWord limit = 1 << n;
        Lut perm;
        perm.reserve(limit);
        if (d <n){
            {
                perm.push_back(0);
                for (int h = 1; h <= d && h <= n; ++h)
                {
                    for (auto x : cpp_vectors_of_hamming_weight(h, n))
                    {
                        perm.push_back(x);
                    }
                }
            }
        }
        BinWord bound=perm.size();
        for (BinWord x = 0; x < limit; ++x)
        {   if (std::popcount(x) > d)
                perm.push_back(x);
        }
        return {std::move(perm),bound};
    }
}

std::vector<cpp_S_box> cpp_all_invariants_up_to_degree(const cpp_S_box S, BinWord d){
    // n=S.get_input_length()
    // debut=time()
    // perm=my_permutation(n,d)
    // inv_perm=perm.inverse()
    // print("Temps passé à calculer perm  et inv_perm",time()-debut)
    // B_S=BinLinearBigBasis([apply_permutation(perm, anf_component(b)) for b in basis_invariants_bis(S)],2**n)
    // bound=sum([binomial(n,t) for t in range(0,d+1)])
    // res=[]
    // for b in B_S.basis_vectors() :
    //     if last_non_zero_index(b) >= bound:
    //         break
    //     res.append(get_sbox(anf_component(apply_permutation(inv_perm,b)))) ## The Mobius Transform is involutive
    // return res
    BinWord n=S.get_input_length();
    auto [perm,bound]= yann_permutation(d,n);
    Lut inv_perm= cpp_inverse(perm);
    cpp_BinLinearBigBasis basis = cpp_BinLinearBigBasis(1 << n);
    for (auto b : cpp_basis_invariants_from_cycles(S)){
        basis.add_to_span(apply_perm(mobius_transform(b,n),perm));
    }
    for (auto b : basis)

}