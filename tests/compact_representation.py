from sage.all import *
from sboxU import *

def are_ea_equivalent_from_vq(f,g):

    # !! TODO !!
    # Add base cases

    ws_f = get_WalshZeroesSpaces(f)
    mappings_f = list(ws_f.get_mappings())


    ws_g = get_WalshZeroesSpaces(g)
    mappings_g = list(ws_g.get_mappings())

    

    quad_index_f = 0
    for i in range(len(mappings_f)):
        f_eq = ccz_equivalent_function(f, mappings_f[i])
        if algebraic_degree(f_eq) == 2:
             quad_index_f = i
             break

    quad_index_g = 0
    for i in range(len(mappings_g)):
        g_eq = ccz_equivalent_function(g, mappings_g[i])
        if algebraic_degree(g_eq) == 2:
             quad_index_g = i
             break
    
    
    # Compute of the ccz-quadratic aut(q)
    q = ccz_equivalent_quadratic_function(f)
    aut_q = automorphisms_from_ortho_derivative(q)
    
    # Compute the EA mapping
    q_prime = ccz_equivalent_quadratic_function(g)
    ea_mappings = ea_mappings_from_ortho_derivative(q_prime,q)
    
    if len(ea_mappings) != 0 : 
        ea_q_q_prime = ea_mappings[0]
    else:
         return(False)

    # To check
    #print(q_prime == ccz_equivalent_function(q,ea_q_q_prime))

    ws_1 =  ws_f.image_by(mappings_f[quad_index_f].inverse().transpose())
    ws_1 =  ws_1.image_by(mappings_f[quad_index_f].inverse().transpose())

    ws_2 =  ws_g.image_by(mappings_g[quad_index_g].inverse().transpose())
    ws_2 =  ws_2.image_by(mappings_g[quad_index_g].inverse().transpose())
    ws_2 =  ws_2.image_by(ea_q_q_prime.transpose())

    #print(ws_1.get_bases()[quad_index_f])
    #print(ws_2.get_bases()[quad_index_g])


    for L in aut_q:
        ws_temp = ws_1.image_by(L.transpose().inverse())
        if ws_2.get_bases()[quad_index_g]  == ws_temp.get_bases()[quad_index_f]:
            return(True)

    return(False)

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





