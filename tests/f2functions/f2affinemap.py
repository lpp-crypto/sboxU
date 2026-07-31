#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *


def main_test():
    with Experiment('Affine Functions of F_2^n'):
        section('Basic functionalities of the F2AffineMap class')
        subsection('Construction')
        # --- { 
        m_mat = Matrix(GF(2), 3, 3, [
            [1, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
        ])
        print(m_mat)
        # --- } 
        # --- { 
        m1 = get_F2AffineMap(m_mat)
        # --- } 
        # --- { 
        m2 = get_F2AffineMap([1, 4, 2])
        # --- } 
        # --- { 
        print("x\tm1(x)\tm2(x)\tmat*x")
        for x in range(0, 2**3):
            print("{}\t{}\t{}\t{}".format(
            x, 
            m1(x),
            m2(x),
            from_bin(m_mat * vector(to_bin(x, 3)))))
        # --- } 
        section('References')
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
