# Code written by Alain Couvreur for
#
#     Canteaut, Anne, Alain Couvreur, and Léo Perrin. "Recovering or
#     testing extended-affine equivalence." IEEE Transactions on
#     Information Theory 68.9 (2022): 6187-6206.
#
# open access version: https://eprint.iacr.org/2021/225


from sage.all import *
from .diff_lin import algebraic_normal_form
from .linear import *


# !SECTION! Some functions on matrices
# ------------------------------------

def right_mult_mat(M):
    # Given a k x n matrix M.
    # returns the (kn) x k^2 matrix representing
    # the right multiplication by M,
    # i.e. the map Mat_{k x k} --> Mat_{k x n}
    #                   X      |-> XM
    k = M.nrows()
    I = identity_matrix(M.base_ring(), k)
    return I.tensor_product(M.transpose(), subdivide = False)


def left_mult_mat(M):
    # Given a k x n matrix M.
    # returns the (kn) x k^2 matrix representing
    # the left multiplication by M,
    # i.e. the map Mat_{n x n} --> Mat_{k x n}
    #                   X      |-> MX
    n = M.ncols()
    I = identity_matrix(M.base_ring(), n)
    return M.tensor_product(I, subdivide = False)


def rank_table(M):
    # Given a k x n matrix M with entries
    # in a polynomial ring in t variables,
    # computes the rank table for any possible
    # entry
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



        

# !SECTION!  Construction of affine equivalent functions
# ------------------------------------------------------


def EAE_function(F, A, a, B, b, C):
    # F a vector function with n inputs and m outputs
    # A an invertible m x m matrix, a an m x 1 matrix
    # B an invertible n x n matrix, b an n x 1 matrix
    # C a full rank m x n matrix
    #
    # returns A*F(Bx + b) + Cx + a           

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





# !SECTION!  Main functions to detect equivalences
# ------------------------------------------------


def solve_system(JF, JG, vv, ww):
    # This function is a single interation of the algorithm
    # Computes (if exist) matrices A, B such that
    # 1. A JG(vi) = JF(w_i) B for vi in vv; wi in ww
    # 2. wi = Bvi for i \geq 2 
    
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
    # Decide the strategy for the enumeration phase
    # which set of values we will guess?
    # The function returns a list of pairs
    # [(r1, m1),...,(rs,ms)] where ri is a rank
    # and mi the number of elements of rank ri
    # we have to guess
    
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
    # returns the remainder of a tuple (a,b,C)
    # from the knowledge of (A,B).
    # The tuple is not unique. We choose the one with b = 0
    m, n = A.ncols(), A.nrows()
    Zm = zero_matrix(GF(2), m, 1)
    b_cand = zero_matrix(GF(2), n, 1)
    G1 = affine_equivalent_function(F, A, Zm, B, b_cand)
    a_cand = (G - G1)([0 for i in range(n)])
    Cpoly = G-G1-a_cand
    C_cand = JacobianMatrix(Cpoly)
    return (a_cand, b_cand, C_cand)




def Alains_get_equivalence(F, G, limit = 10, verbose = False):
    # Function to extract if exist the pair (A,B) of
    # matrices, corresponding to linear parts of the
    # affine equivalence

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
    


# !SECTION! Wrapping Alain's code
# -------------------------------

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
    first_result = Alains_get_equivalence(alainize(f), alainize(g), limit=limit, verbose=verbose)
    return [desAlainize(z) for z in first_result]

