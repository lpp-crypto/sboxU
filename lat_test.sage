from sboxU import *

def scalar_product(alpha, x, p, m):
    sol = 0
    gen = 1
    temp = x
    for i in range(m):
        sol += gen*((alpha*temp)%p)
        temp //= p
        gen *= p
    return sol

def dot_product(a, b, p, m):
    sol = 0
    for i in range(m):
        a_mod_p = a%p
        a = a//p
        b_mod_p = b%p
        b = b//p
        sol += a_mod_p * b_mod_p
    a_mod_p = a%p
    b_mod_p = b%p
    sol += a_mod_p * b_mod_p
    return sol % p


def naive_lat(tab, p, m):
    omega = exp(2*I*pi/p)
    q = p**m
    sol = [[0 for b in range(q)] for a in range(q)]
    for a in range(q):
        for b in range(q):
            temp = 0
            for x in range(q):
                temp += omega ** (dot_product(a,x,p,m) - dot_product(b, tab[x], p, m))
            sol[a][b] = float(abs(temp))
    return sol

if __name__ == '__main__':
    p = 19
    m = 1
    q = p**m
    tab = list(range(q))
    shuffle(tab)

    l1 = lat(tab, p=int(p), n_threads=2)
    l2 = naive_lat(tab, p, m)

    lc1 = lat_column(tab, 17, p=p)
    lc2 = [l2[i][int(17)] for i in range(q)]

    lr1 = lat_row(tab, 17, p=p)
    lr2 = [l2[int(17)][i] for i in range(q)]

    import itertools
    m1 = lat_max(tab, p=p)
    m2 = float(max(l2[i][j] for i,j in itertools.product(range(1,q), range(q))))
    print(m1, m2, sep='\n')

    print(lc1)
    print(lc2)
    print()
    print(lr1)
    print(lr2)

    

    for a in range(q):
        for b in range(q):
            if abs(float(l1[int(a)][int(b)]) - float(l2[a][b])) > 1e-5:
                print('ouch')
            # print(float(l1[a][b] - l2[a][b]))


