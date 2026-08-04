#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *


def main_test():
    with Experiment('Affine Functions of F_2^n'):
        section('Simple Functions')
        section('Basic functionalities of the F2AffineMap class')
        subsection('Construction')
        # --- { 
        m_mat = Matrix(GF(2), 3, 3, [
            [1, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
        ])
        # --- } 
        # --- { 
        m1 = get_F2AffineMap(m_mat)
        # --- } 
        # --- { 
        m2 = get_F2AffineMap([1, 4, 2])
        # --- } 
        # --- { 
        m3 = get_F2AffineMap(get_sbox([0, 1, 4, 5, 2, 3, 6, 7]))
        # --- } 
        # --- { 
        z  = PolynomialRing(GF(2), "z").gen()
        gf = GF(8, modulus=z**3+z+1, name="a")
        a, X  = gf.gen(), PolynomialRing(gf, "X").gen()
        m4 = get_F2AffineMap((a**2+1)*X**4 + (a**2+a+1)*X**2 + (a+1)*X)
        # --- } 
        # --- { 
        print("x\tmat*x\tm1(x)\tm2(x)\tm3(x)\tm4(x)")
        for x in range(0, 2**3):
            y = from_bin(m_mat * vector(to_bin(x, 3)))
            row = "{}\t{}\t".format(x, y)
            for transformation in [m1, m2, m3, m4]:
                y_prime = transformation(x)
                row += "{}\t".format(y_prime)
                if y_prime != y:
                    fail("an F2AffineMap doesn't match the SAGE matrix")
            print(row)
        # --- } 
        subsection('A bigger test: rank distribution')
        # --- { 
        n_max = 10
        parameters = [(7, 10), (10, 7), (10, 10)]
        n_tested = 2**14
        prng = InsecurePRNG(b"seed")
        # --- } 
        # --- { 
        def proba_full_rank(n, m):
            if m > n:
                m, n = n, m
            return float(prod(1 - 2**(-k) for k in range(n - m + 1, n + 1)))
        # --- } 
        # --- { 
        for n, m in parameters:
            counters = [0 for r in range(0, n_max+1)]
            for t in range(0, n_tested):
                L = rand_linear_function(prng, n, m)
                counters[L.rank()] += 1
            row = "({:2d}, {:2d})".format(n, m)
            for c in counters:
                row += "\t{:5.3f}".format(float(c) / n_tested)
            print(row)
        # --- } 
        # --- { 
            observed = float(counters[min(n, m)]) /  n_tested
            expected = proba_full_rank(n, m)
            diff = abs(expected - observed)
            if diff > 0.01:
                fail("mismatch between theory and practice: {}".format(diff))
            
        # --- } 
        section('References')
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
