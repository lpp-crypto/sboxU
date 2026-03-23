
from sboxU.ccz import get_WalshZeroesSpaces,ccz_equivalent_function
from sboxU.apn import ccz_equivalent_quadratic_function, ea_mappings_from_ortho_derivative, automorphisms_from_ortho_derivative
from sboxU.core import algebraic_degree


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
