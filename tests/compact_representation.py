from sage.all import *
from sboxU import *


with Experiment("Testing APN functions-related functions"):

    section("Initialization")
    N = 8
    field = GF(2**8)
    X = PolynomialRing(field, "X").gen()
    g = field.gen()
    
    
    s_boxes  = [ get_sbox(X**3 + X**17 + g**16*(X**18 + X**33) + g**15*X**48),
        get_sbox(X**3 + g**24*X**6 + g**182*X**132 + g**67*X**192),
        get_sbox(X**3 + X**6 + X**68 + X**80 + X**132 + X**160),
        get_sbox(X**3 + X**5 + X**18 + X**40 + X**66),
        get_sbox(X**3 + X**12 + X**40 + X**66 + X**130 ),]


    section("Moebius Transform Check")
    print("The function anf_component computes the Moebius transform, hence its involutive: ")
    a = anf_component(s_boxes[0].coordinate(0))
    b = anf_component(a)
    print(b == s_boxes[0].coordinate(0).lut())
    
    section("Computing The Compact Representation")
    valve = True
    for s in s_boxes:
        qcr = quadratic_compact_representation(s)
        print(qcr)
        lut_from_compact = quadratic_sbox_from_compact_representation(qcr,8,8)
        print("The Compact Representation has correct degree:", end =" ")
        print(algebraic_degree(lut_from_compact) == 2)
        print("The Compact Representation is APN:", end =" ")
        print(is_differential_uniformity_smaller_than(lut_from_compact,2))
        b = are_ea_equivalent_from_vq(s,lut_from_compact)
        print("The Compact Representation is EA-equivalent to the original:", end =" ")
        print(b)
        print()
        if not(b) :
            valve = False
    if valve:
        pprint("[green]All Test Passed")
    else:
        pprint("[red]Error in Conversion")





