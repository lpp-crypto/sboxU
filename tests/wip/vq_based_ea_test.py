# sent on 2025/12/19 on mattermost by Pierre in a private conversation

from sage.all import *
from sboxUv2 import *
from time import *



def are_ea_equivalent_from_vq(f,g):

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
     

n = 6
with APNFunctions('apn6.db') as db:

    with Experiment("testing new EA algorithm"):
        all_apn = db.query_functions({})        
    
    
        section("Generating functions")

        to_test = {}

        subsection("basic, unremarkable functions")
        i = randint(0,716)
        j = randint(0,716)
        s1 = all_apn[i]["sbox"]
        s2 =all_apn[j]["sbox"]
        A = rand_linear_permutation(n)
        B = rand_linear_permutation(n)
        C = rand_linear_function(n, n)
        s3 = Sb(B) * s1 * Sb(A) + Sb(C)
        to_test[ "basic, unremarkable functions" ] = [s1, s2, s3]
        

        subsection("the big automorphism group case")
        s1 = all_apn[0]["sbox"]
        s2 =all_apn[1]["sbox"]
        A = rand_linear_permutation(n)
        B = rand_linear_permutation(n)
        C = rand_linear_function(n, n)
        s3 = Sb(B) * s1 * Sb(A) + Sb(C)

        to_test[ "the big automorphism group case" ] = [s1, s2, s3]
            

        for fs_group in to_test.keys():
            section(fs_group)

            s = to_test[fs_group]
            s1, s2, s3 = s[0], s[1], s[2]

            subsection("Testing a priori in-equivalent functions")
            
            c = Chronograph("using the new algorithm")
            print(are_ea_equivalent_from_vq(s1,s2))
            pprint(c)
            
            c = Chronograph("using *old* algorithm")
            print(are_ea_equivalent(s1,s2))
            pprint(c)
            
            subsection("testing EA-equivalent functions")
            c = Chronograph("using the new algorithm")
            print(are_ea_equivalent_from_vq(s1,s3))
            pprint(c)
            
            c = Chronograph("using *old* algorithm")
            print(are_ea_equivalent(s1,s3))
            pprint(c)
    
