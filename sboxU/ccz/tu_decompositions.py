from sage.all import log
from sboxU.core import is_permutation,identity_F2AffineMap,zero_F2AffineMap,oplus,block_diagonal_F2AffineMap,rank_of_vector_set, get_sbox, get_F2AffineMap,F2AffineMap
from sboxU.ccz import get_WalshZeroesSpaces
from sboxU.algorithms import generating_F2AffineMap_r,generating_F2AffineMap, F2AffineMap_from_range_and_image
from sboxU.display import pprint
from sboxU.statistics import lat


class TUdecomposition:
    """Stores the TU_t decomposition of a function."""
    def __init__(self, A=None, B=None, C=None, T=None, U=None, U_prime=None,m=None):
        if T == None or (U == None and U_prime == None):
            raise Exception("Insufficient information to build a TUdecomposition")
        # setting integer parameters
        self.t = T[0].get_input_length()
        self.n = self.t + int(log(len(T), 2))
        self.mask_t = sum(int(1 << i) for i in range(0, self.t))
        if m is None:
            self.m=self.n
        else :
            self.m=m
        # if U != None :
        #     self.m=U[0].get_output_length()+self.t
        # else :
        #     self.m=U_prime[0].get_output_length()+self.t

        # setting linear components
        if A is None:
            self.A = identity_F2AffineMap(self.n)
        elif isinstance(A, F2AffineMap):
            self.A = A
        else:
            self.A = get_F2AffineMap(A)  

        if B is None:
            self.B = identity_F2AffineMap(self.m)
        elif isinstance(B, F2AffineMap):
            self.B = B
        else:
            self.B = get_F2AffineMap(B)
        if C is None:
            self.C = zero_F2AffineMap(self.n,self.m)
        elif isinstance(C, F2AffineMap):
            self.C = C
        else:
            self.C = get_F2AffineMap(C)
        # setting mini-block ciphers
        self.T = T
        self.T_inv = []
        for T_k in self.T:
            if not is_permutation(T_k):
                raise Exception("ill formed TU-decomposition: T must be a keyed permutation!")
                break
            else:
                self.T_inv.append(get_sbox(T_k.inverse()))
        if U != None:
            self.U = U
            self.U_prime = [[-1 for x2 in range(0, 2**(self.n - self.t))]
                            for x1 in range(0, 2**self.t)]
            for x1 in range(0, 2**self.t):
                for x2 in range(0, 2**(self.n - self.t)):
                    self.U_prime[self.T[x2][x1]][x2] = self.U[x1][x2]
            self.U_prime=[get_sbox(row) for row in self.U_prime]
            # print("U_prime finished")
        elif U_prime != None:
            self.U_prime = U_prime
            self.U = [[-1 for x2 in range(0, 2**(self.n - self.t))]
                      for x1 in range(0, 2**self.t)]
            for x1 in range(0, 2**self.t):
                for x2 in range(0, 2**(self.n - self.t)):
                    self.U[x1][x2] = self.U_prime[self.T[x2][x1]][x2]
            self.U=[get_sbox(row) for row in self.U]
            # print("U finished")
        else:
            raise Exception("at least U or U' must be specified")


    def core(self):
        return TUdecomposition(
            A=None,
            B=None,
            C=None,
            T=self.T,
            U=self.U,
            m=self.m
        )

    
    def twist(self):
        return TUdecomposition(
            A=None,
            B=None,
            C=None,
            T=self.T_inv,
            U=self.U_prime,
            m=self.m
        )


    def insert_before_T(self, alpha):
        """Returns a new TUdecomposition that is functionally
        equivalent, but where the linear permutation `alpha` is
        applied on the t-bit branch on which the T block is a
        permutation *before* it is called.

        """
        new_A = block_diagonal_F2AffineMap(alpha.inverse(),identity_F2AffineMap(self.n-self.t)) *self.A
        alpha_sb=alpha.get_S_box()
        new_T = [row*alpha_sb for row in self.T]
        new_U = [self.U[alpha(y)] for y in range(0, 2**self.t)]
        return TUdecomposition(
            A=new_A,
            B=self.B,
            C=self.C,
            T=new_T,
            U=new_U,
            m=self.m
        )
    
        
    def insert_after_T(self, alpha):
        """Returns a new TUdecomposition that is functionally
        equivalent, but where the linear permutation `alpha` is
        applied on the t-bit branch on which the T block is a
        permutation *after* it is called.

        """
        new_B = self.B* block_diagonal_F2AffineMap(alpha.inverse(), identity_F2AffineMap(self.m-self.t)) 
        alpha_sb=alpha.get_S_box()
        new_T = [alpha_sb*row for row in self.T]
        return TUdecomposition(
            A=self.A,
            B=new_B,
            C=self.C,
            T=new_T,
            U=self.U,
            m=self.m
        )
    


    def insert_before_U(self, alpha):
        """Returns a new TUdecomposition that is functionally
        equivalent, but where the linear permutation `alpha` is
        applied on the (n-t)-bit branch on which the U block is a
        applied, *before* it is called.

        """
        new_A = block_diagonal_F2AffineMap(identity_F2AffineMap(self.t),alpha.inverse() ) * self.A
        alpha_sb=alpha.get_S_box()
        new_T = [self.T[alpha_sb[y]] for y in range(0, 2**(self.n - self.t))]
        new_U = [row *alpha_sb for row in self.U]
        return TUdecomposition(
            A=new_A,
            B=self.B,
            C=self.C,
            T=new_T,
            U=new_U,
            m=self.m
        )
    
        
    def insert_after_U(self, alpha):
        """Returns a new TUdecomposition that is functionally
        equivalent, but where the linear permutation `alpha` is
        applied on the (m-t)-bit branch on which the U block is
        applied, *after* it is called.

        """
        new_B = self.B * block_diagonal_F2AffineMap(identity_F2AffineMap(self.t),alpha.inverse())
        alpha_sb=alpha.get_S_box()
        new_U = [alpha_sb * row for row in self.U]
        return TUdecomposition(
            A=self.A,
            B=new_B,
            C=self.C,
            T=self.T,
            U=new_U,
            m=self.m
        )
    
        

        
    
    def __str__(self):
        result = "t = {:d}\n".format(self.t)
        result += "\nA=[\n"
        result += self.A.__str__()
        result += "\nT=[\n"
        for row in self.T:
            result += "  " + str(row) + "\n"
        result += "]\nU=[\n"
        for row in self.U:
            result += "  " + str(row) + "\n"
        result +=  "]\n"
        result += "]\nU prime=[\n"
        for row in self.U_prime:
            result += "  " + str(row) + "\n"
        result +=  "]\n"
        result += "\nB=[\n"
        result += self.B.__str__()
        result += "\nC = \n"
        result += self.C.__str__()
        return result
    
    
    def get_S_box(self):
        """Returnis the lookup table of the permutation whose TU-decomposition
        is stored.

        """
        result = []
        for x in range(0, 2**self.n):
            y = self.A(x)
            x2, x1 = y >> self.t, y & self.mask_t
            y1 = self.T[x2][x1]
            y2 = self.U[x1][x2]
            y = oplus(y1 | (y2 << (self.t)) , self.C(y))
            y = self.B(y)
            result.append(y)
        return get_sbox(result)
    
    def __eq__(self, value):
        return self.get_S_box()==value.get_S_box()
    

def thickness(basis,n):
    mask = sum(int(1 << i) for i in range(0,   n))
    return rank_of_vector_set([b & mask for b in basis])


    
def tu_decomposition_from_space_basis(s, basis,n,m):
    """Using the knowledge that v is a subspace of Z_s of dimension n, a
    TU-decomposition (as defined in [Perrin17]) of s is performed and
    the corresponding TUdecomposition instance is returned.

    """
    t = thickness(basis, m)
    if t == 0:
        raise Exception("cannot do a TU decomposition when t=0")
    mask_m  = sum(int(1 << i) for i in range(0,   m))
    mask_t  = sum(int(1 << i) for i in range(0,   t))
    mask_n  = sum(int(1 << i) for i in range(0,   n))
    # reordering the basis of the space to have a basis of the thickness space
    reordered_basis = []
    reordered_projected_basis = []
    old_rank = 0
    for b in basis:
        new_reordered_projected_basis = reordered_projected_basis + [b & mask_m]
        new_rank = rank_of_vector_set(new_reordered_projected_basis)
        if new_rank > old_rank:
            reordered_basis.append(b)
            reordered_projected_basis = new_reordered_projected_basis[:]
            old_rank = new_rank
    if len(reordered_basis) != t:
        raise Exception("invalid basis")
    sanitized_basis = reordered_basis[0:t]
    # ensuring that their is no other non-zero vector in the thickness part
    for b in basis:
        if b not in reordered_basis:
            for coeff in range(0, 2**t):
                b_prime = b
                for i in range(0, t):
                    if ((coeff >> i) & 1) == 1:
                        b_prime = oplus(b_prime, reordered_basis[i])
                if (b_prime & mask_m) == 0:
                    sanitized_basis.append(b_prime)
                    break
    # deducing the linear mappings
    basis_A = [b >> m    for b in sanitized_basis[t:n]]
    basis_C = [b >> m   for b in sanitized_basis[0:t]]
    basis_B = [b & mask_m for b in sanitized_basis[0:t]]
    A = generating_F2AffineMap_r(basis_A, n).transpose()
    B = generating_F2AffineMap(basis_B, m).inverse().transpose()
    C = B.inverse()*(F2AffineMap_from_range_and_image(basis_B,basis_C,n,m).transpose())*A.inverse()
    # recovering T and U
    s_prime = []
    for x in range(0, 2**n):
        y = A.inverse()(x)
        y = s[y]
        y = B.inverse()(y)
        s_prime.append(oplus(y, C(x)))
    T = [[0 for x1 in range(0, 2**t)] for x2 in range(0, 2**(n-t))]
    U = [[0 for x2 in range(0, 2**(n-t))] for x1 in range(0, 2**t)]
    for x1 in range(0, 2**t):
        for x2 in range(0, 2**(n-t)):
            x = x1 |(x2 << t)
            y2, y1 = s_prime[x] >> t, s_prime[x] & mask_t
            T[x2][x1] = y1
            U[x1][x2] = y2
    T=[get_sbox(row) for row in T]
    U=[get_sbox(row) for row in U]
    return TUdecomposition(A=A, B=B, C=C, T=T, U=U,m=m)


    
def get_tu_decompositions(s, walsh_zeroes=None):
    sb=get_sbox(s)
    if walsh_zeroes == None:
        walsh_zeroes = get_WalshZeroesSpaces(sb).get_bases()
    result = []
    m=sb.get_output_length()
    n=sb.get_input_length()
    result = []
    for w in walsh_zeroes:
        t = thickness(w, m)
        if t > 0 and t <= m and t<= n:
            d = tu_decomposition_from_space_basis(sb, w,n,m)
            result.append(d)
    return result

def swap_F2AffineMap(t,n,m):
    resultat=[]
    for i in range(t):
        resultat.append(1<<(n+i))
    for i in range(n-t):
        resultat.append(1<<(t+i))
    for i in range(t):
        resultat.append((1<<i))
    for i in range(m-t):
        resultat.append(1<<(t+n+i))
    return get_F2AffineMap(resultat,n+m,n+m)
    
