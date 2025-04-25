import time as time
#from tqdm import *
from sboxU import *


# APN
quad_apn6 = [
[0, 1, 8, 15, 27, 14, 35, 48, 53, 39, 43, 63, 47, 41, 1, 1, 41, 15, 15, 47, 52, 6, 34, 22, 20, 33, 36, 23, 8, 41, 8, 47, 36, 52, 35, 53, 35, 39, 20, 22, 33, 34, 48, 53, 39, 48, 6, 23, 22, 33, 63, 14, 23, 52, 14, 43, 27, 63, 36, 6, 27, 43, 20, 34] ,
[0, 58, 31, 63, 13, 7, 44, 60, 59, 17, 2, 50, 61, 39, 58, 58, 39, 63, 63, 61, 30, 54, 56, 10, 45, 37, 19, 1, 31, 39, 31, 61, 19, 30, 44, 59, 44, 17, 45, 10, 37, 56, 60, 59, 17, 60, 54, 1, 10, 37, 50, 7, 1, 30, 7, 2, 13, 50, 19, 54, 13, 2, 45, 56] ,
[0, 17, 60, 37, 20, 46, 0, 50, 6, 53, 8, 51, 40, 48, 14, 30, 0, 10, 48, 50, 2, 35, 26, 51, 61, 21, 63, 31, 5, 6, 47, 36, 48, 27, 3, 32, 5, 5, 30, 22, 15, 6, 14, 15, 0, 34, 41, 3, 26, 42, 37, 29, 57, 34, 46, 61, 30, 12, 19, 9, 7, 62, 34, 19] ,
[0, 53, 39, 35, 55, 58, 8, 52, 23, 14, 56, 16, 31, 62, 40, 56, 54, 7, 50, 50, 50, 59, 46, 22, 59, 38, 55, 27, 0, 37, 20, 0, 11, 60, 26, 28, 56, 55, 49, 15, 58, 33, 35, 9, 54, 21, 55, 37, 33, 18, 19, 17, 33, 42, 11, 49, 10, 21, 48, 30, 53, 18, 23, 1] ,
[0, 2, 52, 16, 25, 21, 54, 28, 37, 12, 35, 44, 16, 55, 13, 12, 25, 2, 1, 60, 28, 9, 31, 44, 53, 5, 31, 9, 28, 34, 45, 53, 34, 40, 29, 49, 21, 17, 49, 19, 26, 59, 23, 16, 1, 46, 23, 30, 5, 22, 22, 35, 46, 51, 38, 29, 52, 12, 21, 11, 51, 5, 9, 25] ,
[0, 1, 9, 52, 6, 1, 13, 54, 42, 16, 39, 33, 14, 50, 1, 1, 0, 31, 60, 31, 58, 35, 4, 33, 49, 21, 9, 17, 41, 11, 19, 13, 16, 49, 14, 19, 28, 59, 0, 27, 29, 7, 7, 33, 51, 47, 43, 11, 1, 62, 42, 41, 49, 8, 24, 29, 23, 19, 56, 0, 5, 7, 40, 22] ,
[0, 8, 50, 42, 4, 0, 24, 12, 3, 39, 56, 12, 33, 9, 52, 12, 59, 56, 41, 58, 5, 10, 57, 38, 11, 36, 16, 47, 19, 48, 38, 21, 8, 50, 40, 2, 63, 9, 49, 23, 35, 53, 10, 12, 50, 40, 53, 63, 63, 14, 63, 30, 50, 15, 28, 49, 39, 58, 46, 35, 12, 29, 43, 42] ,
[0, 50, 59, 44, 53, 14, 33, 63, 61, 28, 31, 27, 34, 10, 47, 34, 2, 30, 21, 44, 61, 40, 5, 53, 58, 53, 52, 30, 47, 41, 14, 45, 3, 11, 62, 19, 50, 51, 32, 4, 28, 7, 56, 6, 7, 21, 12, 59, 44, 10, 61, 62, 23, 56, 41, 35, 54, 3, 62, 46, 39, 27, 0, 25] ,
[0, 1, 41, 31, 37, 52, 22, 48, 28, 39, 29, 17, 32, 11, 59, 39, 19, 27, 41, 22, 43, 51, 11, 36, 38, 20, 52, 49, 7, 37, 15, 26, 52, 39, 10, 46, 26, 25, 62, 10, 4, 45, 18, 12, 51, 10, 63, 49, 22, 12, 59, 22, 37, 47, 18, 47, 15, 47, 10, 29, 37, 21, 58, 61] ,
[0, 57, 9, 18, 22, 34, 39, 49, 29, 28, 30, 61, 34, 46, 25, 55, 5, 44, 47, 36, 29, 57, 15, 9, 28, 13, 60, 15, 45, 49, 53, 11, 7, 9, 0, 44, 1, 2, 62, 31, 62, 8, 51, 39, 17, 42, 36, 61, 3, 29, 39, 27, 11, 24, 23, 38, 62, 24, 16, 20, 31, 52, 9, 0] ,
[0, 52, 13, 5, 9, 20, 52, 21, 45, 40, 38, 31, 39, 11, 28, 12, 21, 9, 9, 41, 20, 33, 56, 49, 53, 24, 47, 62, 55, 51, 29, 37, 45, 26, 14, 5, 18, 12, 1, 35, 44, 42, 9, 51, 16, 63, 5, 22, 27, 4, 41, 10, 44, 26, 46, 36, 23, 57, 35, 49, 35, 36, 39, 28] ,
[0, 45, 13, 10, 51, 58, 39, 4, 10, 48, 33, 49, 62, 32, 12, 56, 29, 37, 7, 21, 62, 34, 61, 11, 8, 39, 52, 49, 44, 39, 9, 40, 18, 37, 36, 57, 45, 62, 2, 59, 63, 31, 47, 37, 7, 3, 14, 32, 30, 60, 63, 55, 49, 55, 9, 37, 44, 25, 43, 52, 4, 21, 26, 33] ,
[0, 17, 28, 34, 59, 2, 20, 2, 1, 20, 8, 50, 56, 5, 2, 16, 35, 13, 50, 51, 11, 13, 41, 0, 50, 24, 54, 51, 24, 26, 47, 2, 45, 33, 47, 12, 12, 40, 61, 54, 9, 1, 30, 57, 42, 10, 14, 1, 22, 37, 25, 5, 36, 63, 24, 44, 34, 21, 56, 32, 18, 13, 59, 11] ,
]
non_quad_apn6 = [[0, 9, 37, 26, 43, 33, 10, 54, 54, 49, 58, 11, 57, 61, 4, 54, 4, 49, 60, 63, 53, 3, 9, 9, 49, 10, 28, 17, 24, 32, 4, 10, 11, 32, 60, 33, 17, 57, 33, 63, 57, 28, 11, 24, 37, 3, 37, 53, 60, 43, 28, 61, 63, 43, 24, 58, 3, 26, 26, 53, 58, 32, 17, 61],
[0, 9, 33, 43, 26, 37, 10, 54, 54, 49, 58, 11, 57, 61, 4, 54, 4, 49, 60, 63, 53, 3, 9, 9, 49, 10, 32, 24, 17, 28, 4, 10, 11, 33, 60, 17, 33, 57, 32, 63, 57, 37, 11, 37, 24, 3, 28, 53, 60, 24, 43, 61, 63, 28, 43, 58, 3, 17, 32, 53, 58, 26, 26, 61],
[0, 6, 0, 47, 11, 1, 51, 26, 1, 48, 6, 53, 23, 2, 22, 19, 36, 60, 24, 19, 7, 19, 17, 28, 8, 43, 55, 52, 54, 49, 29, 40, 38, 16, 17, 8, 27, 33, 46, 49, 29, 56, 3, 34, 25, 24, 59, 44, 16, 56, 33, 28, 5, 33, 36, 31, 34, 21, 62, 47, 14, 29, 60, 27],
[0, 47, 0, 27, 24, 51, 46, 49, 10, 40, 16, 38, 4, 0, 40, 56, 14, 48, 29, 41, 25, 31, 62, 50, 0, 21, 35, 28, 43, 36, 60, 25, 51, 2, 19, 22, 9, 60, 31, 30, 19, 47, 41, 1, 63, 37, 51, 61, 49, 17, 2, 40, 4, 28, 3, 17, 21, 30, 22, 55, 28, 13, 43, 16],
[0, 51, 54, 57, 37, 46, 31, 40, 37, 14, 38, 63, 43, 54, 0, 47, 61, 44, 38, 5, 14, 39, 29, 6, 8, 1, 58, 15, 10, 53, 20, 23, 21, 40, 3, 2, 52, 49, 46, 23, 54, 19, 21, 2, 60, 47, 55, 22, 50, 45, 9, 36, 5, 34, 54, 35, 1, 6, 19, 40, 7, 54, 57, 52]]




##############

# Computes the lut of f+u*eps with f a lut of int, u an int and eps of vector of boolean
def switch(f,u,eps,n):
    res = []
    for i in range(2**n):
        if int(eps[i]) ==1 :
            res.append(oplus(f[i],u))
        else :
            res.append(f[i])
    return(res)    
    
def diff(f,a,x):
    return(oplus(f[oplus(a,x)],f[x]))


# Computes naively all the switching criteria for lut
def naive_switching_criterion(lut,n):
    BF = GF(2)
    constraint_list = [[] for i in range(2**n)]
    # For each u, we check naively if T(x+a)+T(x)+T(y)+T(y+a) = u
    for a in range(2**n) :
        for x in range(2**n) :
            for y in range(2**n) :
                u = oplus(diff(lut,a,x),diff(lut,a,y))
                temp = [BF(0)]*(2**n)
                temp[oplus(x,a)] = BF(1)
                temp[x] = BF(1)
                temp[oplus(y,a)] = BF(1)
                temp[y] = BF(1)
                constraint_list[u].append(temp)
                
    return(constraint_list)


## Affine Basis ##

# Computes the basis of affine functions as vectors of 2**n elements
def compute_affine_basis(n):
    BF = GF(2)
    basis = []
    # Constant
    basis.append([1]*(2**n))
    # The rest
    for i in range(n):
        temp = [ BF(1)*(j>>i)%2 for j in range(2**n)]
        basis.append(temp)
    return(basis)


## Removing Affine Equations

# Computes the complement of span(affine_basis) in the vector space span(switching_criterion)
def remove_affine_equations(switching_criterion_u,affine_basis):

    ker_crit_basis = list(matrix(switching_criterion_u).right_kernel().basis())
    V = VectorSpace(GF(2),64)
    W = V.subspace(ker_crit_basis)
    Wtemp = V.subspace(affine_basis)
    non_affine_basis = []
    # Draws a vector at random until the dimension increases (at worst it is proba 1/2)
    while Wtemp.dimension() < W.dimension() :
            w = W.random_element()
            if (Wtemp+span([w])).dimension() > Wtemp.dimension():
                Wtemp  = Wtemp+span([w])
                non_affine_basis.append(w)

    return(V.subspace(non_affine_basis))




## APN Switching Class


# Computes the apn functions in the switching class of lut
def apn_non_affine_switching_class(lut):
    n = int(log(len(lut), 2))
    switch_list = []
    affine_basis = compute_affine_basis(n)
    switching_criterion = naive_switching_criterion(lut,n)
    for u in range(2**n):
        Veps = remove_affine_equations(switching_criterion[u],affine_basis)
        for eps in Veps:
            switch_list.append(switch(lut,u,eps,n))
    return(switch_list)

# Returns True if there si non quadratic elements in the switching class of lut
def is_there_non_quad_apn_switch(lut,n):
    switch_list = []
    affine_basis = compute_affine_basis(n)
    switching_criterion = naive_switching_criterion(lut,n)
    for u in range(2**n):
        Veps = remove_affine_equations(switching_criterion[u],affine_basis)
        for eps in Veps.basis():
            ANF = algebraic_normal_form([int(b) for b in eps])
            if ANF[0].degree() > 2 :
                return(True)
    return(False)



## Example
"""
def max_matrix(M):
    c = 0
    for i in M[1:]:
        for j in i :
            if j > c :
                c = j
    return(c)

t = time.time()
switch_list  = apn_non_affine_switching_class(quad_apn6[2],6)
for i in switch_list:
    if max_matrix(ddt(i)) > 2 :
        print('oops')
print(time.time() - t)
"""





