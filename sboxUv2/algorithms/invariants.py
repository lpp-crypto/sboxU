from sage.all import Permutation, Combinations, binomial
from sboxUv2.core.sbox import Sb
from sboxUv2.core.f2functions import to_bin, from_bin
from sboxUv2.algorithms import BinLinearBigBasis
from sboxUv2.core.anf import anf_component

def basis_invariants(S):
    """
    Computes a basis of invariant Boolean functions
    associated to the permutation S (TLS19).
    """
    all_even = True
    R = []

    cycles = Permutation([x + 1 for x in S.lut()]).to_cycles()
    n=S.input_space_size()

    for cycle in cycles:
        if len(cycle) % 2 == 1:
            all_even = False

    B = [0] * n
    for cycle in cycles:
        for l in cycle:
            B[l - 1] = 1
        R.append(Sb(B))
        B = [0] * n

    if all_even:
        B = [0] * n
        for cycle in cycles:
            value = 0
            for i in range(len(cycle)):
                B[cycle[i] - 1] = value
                value = (value + 1) % 2
        R.append(Sb(B))

    return R

def int_of_hamming_weigth(h,n):
    res=[]
    for places in Combinations([i for i in range(n)], h):
        x=[0 for _ in range(n)]
        for i in places:
            x[i]=1
        res.append(from_bin(x))
    return res

def my_permutation(n,d):
    perm=[]
    to_assign=[i for i in range(2**n)]
    for h in range(d+1):
        for x in int_of_hamming_weigth(h,n):
            perm.append(x)
            to_assign.remove(x)
    perm+=to_assign
    return Sb(perm)

def last_non_zero_index(v):
    for i in range(len(v)-1,0,-1):
        if v[i]==1:
            return i
    return 0


def apply_permutation(perm,vect):
    return [vect[perm[i]] for i in range(len(vect))]


def all_invariants_up_to_degree(S,d):
    n=S.get_input_length()
    perm=my_permutation(n,d)
    inv_perm=perm.inverse()
    B_S=BinLinearBigBasis([apply_permutation(perm, anf_component(b)) for b in basis_invariants(S)],2**n)
    bound=sum([binomial(n,t) for t in range(0,d+1)])
    res=[]
    for b in B_S.basis_vectors() :
        if last_non_zero_index(b) >= bound:
            break
        res.append(Sb(anf_component(apply_permutation(inv_perm,b)))) ## The Mobius Transform is involutive
    return res