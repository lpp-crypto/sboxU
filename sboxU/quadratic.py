#!/usr/bin/sage
# Time-stamp: <2024-01-05 14:12:51 lperrin>


from .diff_lin import *
from .linear import *

from collections import defaultdict


# !SECTION! QIC
# =============

def get_empty_QIC(N):
    return [[[0 for k in range(0, N)]
             for j in range(0, N)]
            for i in range(0, N)]


def get_QIC(s):
    """Returns the QIC of s, which is essentially the
    linear part of its Jacobian except that each entry is represented by a
    vector rather than a linearized polynomial.

    """
    N = int(log(len(s), 2))
    e = [int(1 << i) for i in range(0, N)]
    Q = get_empty_QIC(N)
    for i in range(0, N):
        s_i = [scal_prod(e[i], s[x]) for x in range(0, len(s))]
        for j in range(0, N):
            for k in range(0, j):
                Q[i][j][k] = oplus(
                    oplus(s_i[oplus(e[j], e[k])], s_i[e[k]]), # \Delta_{e_j}S_i(e_k)
                    oplus(s_i[oplus(e[j], 0)], s_i[0]),       # \Delta_{e_j}S_i(0)
                )
    return Q


def from_QIC(Q):
    """Returns the LUT of a function with QIC `Q`. """
    N = len(Q)
    e = [int(1 << i) for i in range(0, N)]
    result = [0]
    for x in range(1, 2**N):
        y = 0
        for i in range(0, N):
            y_i = 0
            for j in range(0, N):
                if scal_prod(x, e[j]) == 1:
                    for k in range(0, j):
                        y_i = oplus(y_i, scal_prod(x, e[k]) * Q[i][j][k])
            y += int(y_i << i)
        result.append(y)
    return result


# !SECTION! Ortho-integration
# ===========================


def ortho_integration(o):
    """Given the lookup `o` of a function, returns (if it exists) the
    lookup table of quadratic APN function that has `o` as its
    ortho-derivative. If it doesn't exist, returns `None`.

    """
    N = int(log(len(o), 2))
    e = [int(1 << i) for i in range(0, N)]
    if o[0] != 0:
        raise Exception("ortho-derivative must map 0 to 0!")
    # building the system that the QIC must be solution of
    row_length = int(N**3)
    mat = []
    for a in range(1, 2**N):
        # ensuring that, for all a, scal_prod(pi_F(a), J(a)(e_k)) = 0
        for k in range(0, N):
            row = [0 for l in range(0, row_length)]
            for i in range(0, N):
                o_i = scal_prod(e[i], o[a])
                if o_i > 0:
                    for j in range(0, N):
                        if scal_prod(e[j], a) > 0:
                            row[i*N**2 + j*N + k] = 1
                            row[i*N**2 + k*N + j] = 1
            mat.append(row)
    # adding contraints to enforce that Q^i_{j,k} = 0 if j <= k
    for i in range(0, N):
        for j in range(0, N):
            for k in range(j, N):
                row = [1*(l == i*N**2 + j*N + k) for l in range(0, row_length)]
                mat.append(row)
    # solving the system
    mat = Matrix(GF(2), len(mat), row_length, mat)
    sol = mat.right_kernel()
    result = []
    # for each solution, rebuild the QIC, and then compute the
    # corresponding function
    for assignment in sol.basis():
        assignment = list(assignment)
        Q = get_empty_QIC(N)
        # removing the all-zero solution
        valid_solution = False
        for x in assignment:
            if x != 0:
                valid_solution = True
                break
        # recovering the QIC and then the function
        if valid_solution:
            for i in range(0, N):
                for j in range(0, N):
                    for k in range(0, N):
                        Q[i][j][k] = assignment[i*N**2 + j*N + k]
            result.append(from_QIC(Q))
    return result


def test_ortho_integration():
    six_bit_apns = [
        [0, 1, 8, 15, 27, 14, 35, 48, 53, 39, 43, 63, 47, 41, 1, 1, 41, 15, 15, 47, 52, 6, 34, 22, 20, 33, 36, 23, 8, 41, 8, 47, 36, 52, 35, 53, 35, 39, 20, 22, 33, 34, 48, 53, 39, 48, 6, 23, 22, 33, 63, 14, 23, 52, 14, 43, 27, 63, 36, 6, 27, 43, 20, 34],
        [0, 58, 31, 63, 13, 7, 44, 60, 59, 17, 2, 50, 61, 39, 58, 58, 39, 63, 63, 61, 30, 54, 56, 10, 45, 37, 19, 1, 31, 39, 31, 61, 19, 30, 44, 59, 44, 17, 45, 10, 37, 56, 60, 59, 17, 60, 54, 1, 10, 37, 50, 7, 1, 30, 7, 2, 13, 50, 19, 54, 13, 2, 45, 56],
        [0, 17, 60, 37, 20, 46, 0, 50, 6, 53, 8, 51, 40, 48, 14, 30, 0, 10, 48, 50, 2, 35, 26, 51, 61, 21, 63, 31, 5, 6, 47, 36, 48, 27, 3, 32, 5, 5, 30, 22, 15, 6, 14, 15, 0, 34, 41, 3, 26, 42, 37, 29, 57, 34, 46, 61, 30, 12, 19, 9, 7, 62, 34, 19],
        [0, 53, 39, 35, 55, 58, 8, 52, 23, 14, 56, 16, 31, 62, 40, 56, 54, 7, 50, 50, 50, 59, 46, 22, 59, 38, 55, 27, 0, 37, 20, 0, 11, 60, 26, 28, 56, 55, 49, 15, 58, 33, 35, 9, 54, 21, 55, 37, 33, 18, 19, 17, 33, 42, 11, 49, 10, 21, 48, 30, 53, 18, 23, 1],
        [0, 2, 52, 16, 25, 21, 54, 28, 37, 12, 35, 44, 16, 55, 13, 12, 25, 2, 1, 60, 28, 9, 31, 44, 53, 5, 31, 9, 28, 34, 45, 53, 34, 40, 29, 49, 21, 17, 49, 19, 26, 59, 23, 16, 1, 46, 23, 30, 5, 22, 22, 35, 46, 51, 38, 29, 52, 12, 21, 11, 51, 5, 9, 25],
        [0, 1, 9, 52, 6, 1, 13, 54, 42, 16, 39, 33, 14, 50, 1, 1, 0, 31, 60, 31, 58, 35, 4, 33, 49, 21, 9, 17, 41, 11, 19, 13, 16, 49, 14, 19, 28, 59, 0, 27, 29, 7, 7, 33, 51, 47, 43, 11, 1, 62, 42, 41, 49, 8, 24, 29, 23, 19, 56, 0, 5, 7, 40, 22],
        [0, 8, 50, 42, 4, 0, 24, 12, 3, 39, 56, 12, 33, 9, 52, 12, 59, 56, 41, 58, 5, 10, 57, 38, 11, 36, 16, 47, 19, 48, 38, 21, 8, 50, 40, 2, 63, 9, 49, 23, 35, 53, 10, 12, 50, 40, 53, 63, 63, 14, 63, 30, 50, 15, 28, 49, 39, 58, 46, 35, 12, 29, 43, 42],
        [0, 50, 59, 44, 53, 14, 33, 63, 61, 28, 31, 27, 34, 10, 47, 34, 2, 30, 21, 44, 61, 40, 5, 53, 58, 53, 52, 30, 47, 41, 14, 45, 3, 11, 62, 19, 50, 51, 32, 4, 28, 7, 56, 6, 7, 21, 12, 59, 44, 10, 61, 62, 23, 56, 41, 35, 54, 3, 62, 46, 39, 27, 0, 25],
        [0, 1, 41, 31, 37, 52, 22, 48, 28, 39, 29, 17, 32, 11, 59, 39, 19, 27, 41, 22, 43, 51, 11, 36, 38, 20, 52, 49, 7, 37, 15, 26, 52, 39, 10, 46, 26, 25, 62, 10, 4, 45, 18, 12, 51, 10, 63, 49, 22, 12, 59, 22, 37, 47, 18, 47, 15, 47, 10, 29, 37, 21, 58, 61],
        [0, 57, 9, 18, 22, 34, 39, 49, 29, 28, 30, 61, 34, 46, 25, 55, 5, 44, 47, 36, 29, 57, 15, 9, 28, 13, 60, 15, 45, 49, 53, 11, 7, 9, 0, 44, 1, 2, 62, 31, 62, 8, 51, 39, 17, 42, 36, 61, 3, 29, 39, 27, 11, 24, 23, 38, 62, 24, 16, 20, 31, 52, 9, 0],
        [0, 52, 13, 5, 9, 20, 52, 21, 45, 40, 38, 31, 39, 11, 28, 12, 21, 9, 9, 41, 20, 33, 56, 49, 53, 24, 47, 62, 55, 51, 29, 37, 45, 26, 14, 5, 18, 12, 1, 35, 44, 42, 9, 51, 16, 63, 5, 22, 27, 4, 41, 10, 44, 26, 46, 36, 23, 57, 35, 49, 35, 36, 39, 28],
        [0, 45, 13, 10, 51, 58, 39, 4, 10, 48, 33, 49, 62, 32, 12, 56, 29, 37, 7, 21, 62, 34, 61, 11, 8, 39, 52, 49, 44, 39, 9, 40, 18, 37, 36, 57, 45, 62, 2, 59, 63, 31, 47, 37, 7, 3, 14, 32, 30, 60, 63, 55, 49, 55, 9, 37, 44, 25, 43, 52, 4, 21, 26, 33],
[0, 17, 28, 34, 59, 2, 20, 2, 1, 20, 8, 50, 56, 5, 2, 16, 35, 13, 50, 51, 11, 13, 41, 0, 50, 24, 54, 51, 24, 26, 47, 2, 45, 33, 47, 12, 12, 40, 61, 54, 9, 1, 30, 57, 42, 10, 14, 1, 22, 37, 25, 5, 36, 63, 24, 44, 34, 21, 56, 32, 18, 13, 59, 11]]
    for s in six_bit_apns:
        print("original : ", pretty_spectrum(differential_spectrum(s)))
        o = ortho_derivative(s)
        for p in ortho_integration(o):
            print("\n")
            print(p)
            print(pretty_spectrum(differential_spectrum(p)))
            print("  {}\n  {} (ref)".format(
                pretty_spectrum(differential_spectrum(ortho_derivative(p))),
                pretty_spectrum(differential_spectrum(o))
            ))
            if len(ea_equivalence_quadratic(p, s)) == 0:
                print("[ERROR] incorrect ortho-integral")
            else:
                print("[SUCCESS]")
                   




# !SECTION! Some functions on matrices
# ------------------------------------

def right_mult_mat(M):
    """Given a k x n matrix M.  returns the (kn) x k^2 matrix
     representing the right multiplication by M, i.e. the map

    Mat_{k x k} --> Mat_{k x n}
         X      |-> XM

    """
    k = M.nrows()
    I = identity_matrix(M.base_ring(), k)
    return I.tensor_product(M.transpose(), subdivide = False)


def left_mult_mat(M):
    """Given a k x n matrix M. returns the (kn) x k^2 matrix
    representing the left multiplication by M, i.e. the map
    
    Mat_{n x n} --> Mat_{k x n}
         X      |-> MX

    """
    n = M.ncols()
    I = identity_matrix(M.base_ring(), n)
    return M.tensor_product(I, subdivide = False)


def rank_table(M):
    """Given a k x n matrix M with entries in a polynomial ring in t
    variables, computes the rank table for any possible entry.

    """
    assert(M.nrows() != 0 and M.ncols() != 0), "[rank_table] : input matrix should be non empty"
            
    max_rank = min(M.nrows(), M.ncols())
    rank_table = [[] for i in range(max_rank + 1)]
    R = parent(M[0,0])
    
    for v in VectorSpace(GF(2), R.ngens()):
        r = M(v.list()).rank()
        rank_table[r].append(v)

    return rank_table


def rank_distribution(rt):
    return [len(l) for l in rt]


def JacobianMatrix(F):
    assert(F.ncols() == 1), "[JacobianMatrix] : The input function should be an m x 1 matrix"

    m = F.nrows()
    assert(m != 0), "[JacobianMatrix] : The input function should be non empty"

    R = F.base_ring()
    n = R.ngens()
    rows = [[F[j,0].derivative(xi) for xi in R.gens()] for j in range(m)]
    return Matrix(R, m, n, rows)

        

# !SUBSECTION!  Construction of affine equivalent functions


def EAE_function(F, A, a, B, b, C):
    """Considering that its input are
    `F` a vector function with n inputs and m outputs,
    `A` an invertible m x m matrix, a an m x 1 matrix,
    `B` an invertible n x n matrix, b an n x 1 matrix,
    `C` a full rank m x n matrix,
    
    returns A*F(Bx + b) + Cx + a           

    """
    ## 1. Some tests
    assert(F.ncols() == 1), "[EAE_function] : input F should be a vector function"

    m = F.nrows()
    assert(m != 0), "[EAE_function] : input F should be a non empty function"
    
    assert(A.nrows() == A.ncols()), "[EAE_function] : matrix A should be square"
    
    assert(B.nrows() == B.ncols()), "[EAE_function] : matrix B should be square"

    assert(C.nrows() == A.ncols()), "[EAE_function] : matrix C has not the good number of rows"
    assert(C.ncols() == B.nrows()), "[EAE_function] : matrix C has not the good number of columns"

    
    R = F.base_ring()
    n = R.ngens()

    assert(A.nrows() == m), "[EAE_function] : matrix A's size should be the output size of F"
        
    assert(B.nrows() == n), "[EAE_function] : matrix B's size should be number of variables of F"
    
    assert(C.ncols() == n), "[EAE_function] : matrix C's size should be number of variables of F"
        
    assert(a.nrows() == m and a.ncols() == 1), "[EAE_function] : a should be an m x 1 matrix"
        
    assert(b.nrows() == n and b.ncols() == 1), "[EAE_function] : b should be an n x 1 matrix"
        
    
    ## 2. Construction of the equivalent function
    var = Matrix(R, n, 1, R.gens())
    variable_change = B * var + b
            
    # Construction of the new function
    G0 = F(variable_change.list())
    
    # Reduce modulo the field equations
    I = R.ideal([x**2 - x for x in R.gens()])
    for i in range(m):
        G0[i,0] = G0[i,0].mod(I)
            
    return A * G0 + C * var + a


def affine_equivalent_function(F, A, a, B, b):
    C = zero_matrix(GF(2), A.nrows(), B.nrows())
    return EAE_function(F, A, a, B, b, C)



# !SUBSECTION!  Main functions to detect equivalences


def solve_system(JF, JG, vv, ww):
    """This function is a single interation of the algorithm Computes
    (if exist) matrices A, B such that
    
    1. A JG(vi) = JF(w_i) B for vi in vv; wi in ww
    2. wi = Bvi for i \\geq 2
    
    """
    ## 1. Conditions
    head_str= "[system_of_matrices] : "
    assert(JG.nrows() == JF.nrows()), head_str + "JF and JG should have the same number of rows"
    assert(JF.ncols() == JG.ncols()), head_str + "JF and JG should have the same numer of columns"
    assert(len(vv) == len(ww)), head_str + "vv and ww should have the same length"
    assert(len(vv) != 0), head_str + "vv and ww should be nonempty"
        
    R = JF.base_ring()
    assert(JG.base_ring() == R), head_str + "JF and JG should have the same number of variables"

    m = JF.nrows()
    n = R.ngens()
    for v in vv:
        assert(len(v) == n), head_str + "entries of vv are lists of length n"
        
    for w in ww:
        assert(len(w) == n), head_str + "entries of ww are lists of length n"
    

    ## 2. Construction of the matrices we need            
    LJF = [left_mult_mat(JF(w.list())) for w in ww]
    RJG = [right_mult_mat(JG(v.list())) for v in vv]
    Rv = [right_mult_mat(Matrix(GF(2), n, 1, v.list()))
          for v in vv]

    List = [[RJG[i], -LJF[i]] for i in range(len(LJF))]
    List += [[zero_matrix(GF(2),n, m**2), r] for r in Rv]
    M = block_matrix(List)
    List2 = [0 for i in range(m * n * len(vv))]
    for w in ww:
          List2 += w.list()
          
    V = matrix(GF(2), (m + 1) * n * len(vv), 1, List2)

    ## 3. Solve the system
    assert(M.ncols() == m**2 + n**2), "[Extract_AB] : M should have n^2 + m^2 columns"
    try:
        T = M.solve_right(V)
    except ValueError:
        # print( "No solution" )
        return None

    return (M.right_kernel(), T, M, V)## remove M, V when done




def AB_from_vector(vect, m, n):
    caution = "[AB_from_vector] : first entry should have length m^2+n^2"
    assert(len(vect) == m*m + n*n), caution
    Am1 = matrix(GF(2), m, m, vect[ : m*m])
    B = matrix(GF(2), n, n, vect[m*m : ])
    return (Am1, B)
    


def strategy(rank_tab, rank_dist, m, n):
    """Decide the strategy for the enumeration phase which set of
    values we will guess?  The function returns a list of pairs [(r1,
    m1),...,(rs,ms)] where ri is a rank and mi the number of elements
    of rank ri we have to guess

    """
    ## 1. Extracts a sorted distribution of interest
    distrib = rank_dist[:]
    distrib.remove(1)
    dist_set = set(distrib)
    dist_set.remove(0)
    distrib = list(dist_set)
    distrib.sort()

    ## 2. Defines some variables
    mx = max(m, n)
    n_equations = 0
    n_variables = m**2 + n**2
    reference_vectors = []
    strat = []

    ## 3. Go!
    max_rank = min(m, n)
    for d in distrib:
        rank = max_rank

        ## 3.1. Finds the right most entry of cardinality d in the
        ## rank table
        while rank > 0 and rank_dist[rank] != d:
            rank -= 1

        ## 3.2. Fills in the output list    
        included_for_this_rank = 0    
        for i in range(d):
            M = matrix(reference_vectors + [rank_tab[rank][i]])
            if M.rank() == len(reference_vectors) + 1:
                reference_vectors.append(rank_tab[rank][i])
                n_equations += rank * mx + n
                included_for_this_rank += 1
                if n_equations >= n_variables:
                    strat.append((rank, included_for_this_rank))
                    return (reference_vectors, strat)

        strat.append((rank, included_for_this_rank))        
    return None



def check_candidate(JF, JG, A, B):
    n = A.nrows()
    for v in basis(VectorSpace(GF(2), n)):
        if JG(v.list()) != A * JF((B * v).list()) * B:
            return False

    return True



def full_tuple(F,G,A,B):
    """Returns the remainder of a tuple (a,b,C) from the knowledge of
    (A,B).  The tuple is not unique. We choose the one with b = 0.

    """
    m, n = A.ncols(), A.nrows()
    Zm = zero_matrix(GF(2), m, 1)
    b_cand = zero_matrix(GF(2), n, 1)
    G1 = affine_equivalent_function(F, A, Zm, B, b_cand)
    a_cand = (G - G1)([0 for i in range(n)])
    Cpoly = G-G1-a_cand
    C_cand = JacobianMatrix(Cpoly)
    return (a_cand, b_cand, C_cand)




def Alains_get_equivalence(F, G, limit = 10, verbose = False):
    """Function to extract if exist the pair (A,B) of matrices,
    corresponding to linear parts of the affine equivalence.

    """

    ## 1. Jacobian matrices and verifications
    JF1 = JacobianMatrix(F)
    JG1 = JacobianMatrix(G)
    n = JG1.ncols()
    m = JG1.nrows()

    caution = "[get_extended_affine_equivalence] : "
    assert(JF1.ncols() == n), caution + "F, G should have the same number of variables"
    assert(JF1.nrows() == m), caution + "F, G should have the same number of outputs"

    ## 2. The linear part of the jacobian matrices
    JF = JF1 - JF1([0 for i in range(n)])
    JG = JG1 - JG1([0 for i in range(n)])
    rkt1 = rank_table(JF)
    rkt2 = rank_table(JG)
    rank_dist1 = rank_distribution(rkt1)
    rank_dist2 = rank_distribution(rkt2)
    if verbose:
        print( "\n------------" )
        print( "Rank distributions : " )
        print( rank_dist1 )
        print( rank_dist2 )
        print( "------------\n" )

    ## 3. Compares the rank distributions 
    if rank_dist1 != rank_dist2:
        print( "F and G have not the same Jacobian rank distribution" )
        print( "They are not equivalent!\n" )
        return None

    ## 4. Defines the strategy of enumeration
    STRAT = strategy(rkt2, rank_dist2, m, n)
    reference_vectors = STRAT[0]
    strat = STRAT[1]

    iter_list = []
    for pair in strat:
        l = len(rkt1[pair[0]])
        iter_list.append(Arrangements(range(l), pair[1]))

    
    ## 5. Enumeration
    ww = [VectorSpace(GF(2), n)(0) for i in range(len(reference_vectors))]
    counter = 0
    
    for index_list in cartesian_product(iter_list):
        counter += 1
        i = 0
        ## 5.1. Construction of the tuple ww of candidates
        ## to be the B*v_i's
        for k in range(len(strat)):
            for index in index_list[k]:
                ww[i] = rkt1[strat[k][0]][index]
                i += 1

        if verbose:        
            print( "-----" )
            print( str(ww) )

        if matrix(ww).rank() == len(ww):
            Sol = solve_system(JF, JG, reference_vectors, ww)
            if Sol == None:
                if verbose:
                    print( "No solution for this tuple" )
                continue

            if dimension(Sol[0]) > limit:
                if verbose:
                    print( "One candidate with too many solutions" )
                ## To do -- traiter ce cas proprement -> Appel supplementaire
                continue

            x0 = vector(Sol[1])
            if verbose:
                print( "Get in the kernel of dimension : " + str(dimension(Sol[0])) )

            for x in Sol[0]:
                AB =  AB_from_vector(x0 + x, m, n)
                if AB[0].is_invertible():
                    A_cand = AB[0]**(-1)
                    B_cand = AB[1]
                    if check_candidate(JF, JG, A_cand, B_cand):
                        abC = full_tuple(F, G, A_cand, B_cand)
                        a_cand = abC[0]
                        b_cand = abC[1]
                        C_cand = abC[2]
                        G1 = EAE_function(F, A_cand, a_cand, B_cand, b_cand, C_cand)
                        assert(G1 == G), caution + "EAE equivalence not found"
                        if verbose:
                            print( "Equivalence found after " + str(counter) + " tries!" )
                            print( "\n Rank tables : " )
                            print( rank_dist1 )
                            print( "\n" )
                        return (A_cand, a_cand, B_cand, b_cand, C_cand)
    return []
    


# !SUBSECTION! Wrapping Alain's code

def alainize(s):
    """The objects manipulated in Alain's code are differ from the
    rest of sboxU. This function acts as an interface between both
    representations.

    """
    anf = algebraic_normal_form(s)
    n = len(anf)
    R_Alain = PolynomialRing(GF(2), "x", n)
    R_boolean = anf[0].parent()
    correspondance = {R_boolean.gens()[i] : R_Alain.gens()[i] for i in range(0, n)}
    Alain_style_polynomials = []
    for a in anf:
        p = R_Alain.zero()
        for monomial in a:
            new_monomial = R_Alain.one()
            for v in monomial:
                new_monomial = new_monomial * correspondance[v]
            p += new_monomial
        Alain_style_polynomials.append(p)
    return Matrix(len(Alain_style_polynomials), 1, Alain_style_polynomials)


def desAlainize(v):
    """The objects manipulated in Alain's code are differ from the
    rest of sboxU. This function acts as an interface between both
    representations.

    """
    if v.ncols() == 1:
        # case of a vector
        return frombin(list([v[i][0] for i in reversed(range(0, v.nrows()))]))
    else:
        # case of a Matrix
        return Matrix(GF(2), v.nrows(), v.ncols(), [
            [v[v.nrows()-1-i][v.ncols()-1-j] for j in range(0, v.nrows())]
            for i in range(0, v.nrows())
        ])


def ea_equivalence_quadratic(f, g, limit=10, verbose=False):
    """Tests the extended-affine equivalence of the functions with LUT
    `f` and `g` assuming that they are quadratic. Based on the results
    presented in

    Canteaut, A., Couvreur, A., Perrin, L.: Recovering or testing
    extended-affine equiv- alence. IEEE Trans. Inform. Theory 68(9),
    6187–6206 (2022). https://doi.org/10.1109/TIT.2022.3166692

    """    
    first_result = Alains_get_equivalence(alainize(f), alainize(g), limit=limit, verbose=verbose)
    return [desAlainize(z) for z in first_result]
