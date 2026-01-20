from sage.all import *
from sboxUv2 import *
# Utils for display #
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import *

##################################################################################################
# This script computes the functions at distance 1 in the CCZ-SN graph for all known  6-bits APN #
##################################################################################################



def spectrify(sn_f):
    spec = {}
    for u in range(len(sn_f)):
        sn_fu_size = log(len(sn_f[u])+1,2)
        #sn_fu_size = len(sn_f[u])
        if sn_fu_size in spec.keys():
            spec[sn_fu_size] += 1
        else:
            spec[sn_fu_size] = 1
    return(spec)


def at_distance_1_6(f,n,db,known_ccz_id,known_id):

    # Computing Neighborhood 
    sn_f = non_trivial_sn(f["sbox"],200,20)

    # We store the ccz_id and the id of the functions we already went through 
    known_ccz_id.update(set([f["ccz_id"]]))
    # We add all the ea id
    l = db.query_functions({'ccz_id':f["ccz_id"]})
    for ll in l:
        known_id.add(ll["id"])


    #####################
    # Specific to n = 6 #
    #####################

    non_ccz_quad_id = set([t["id"] for t in db.query_functions({'ccz_id': 13})])

    # Db id of the only quadratic functions
    quad_id = set([fun["id"] for fun in db.query_functions({"degree":2})])

    # For each function in the Neighborhood, we check the EA representatives of its CCZ class and check where we land
    for u in range(len(sn_f)):

        for i in range(len(sn_f[u])):
            s = sn_f[u][i]
            mug_s = apn_ea_mugshot(s)
            same_mug_s = db.query_functions({'mugshot': mug_s})
            if same_mug_s == []:
                print("NEW FUNCTION, NEW MUGSHOT")
                print(s)
            new_valve = True
            same_mug_s_ccz_id = set([t["ccz_id"] for t in same_mug_s])
            same_mug_s_id = set([t["id"] for t in same_mug_s])

            # If we have all the possible CCZ classes, no need to continue
            check_set = known_ccz_id.intersection(same_mug_s_ccz_id)
            if not(len(same_mug_s_ccz_id.difference(check_set))== 0):

                # If s is quadratic we can check equivalence against only quadratic functions
                if algebraic_degree(s) == 2 :

                    good_set = quad_id.intersection(same_mug_s_id)
                    # We then test the EL-equivalence of ortho-derivatives
                    for j in good_set :

                        pi_s = ortho_derivative(s)
                        pi_j = ortho_derivative(db[j]["sbox"].lut())
                        if are_linear_equivalent(pi_j,pi_s):

                            known_ccz_id.add(db[j]["ccz_id"])
                            known_id.update(set([fun["id"] for fun in db.query_functions({"ccz_id" : db[j]["ccz_id"]})]))
                            new_valve = False
                            break

                else :
                    s_quad = ccz_equivalent_quadratic_function(s)
                    # If s is CCZ_quadratic, we check the EL of its ortho-derivative of its quadratic representative to a matching function
                    if len(s_quad) != 0 :

                        # We check EL-equivalence of functions that look like s_quad
                        mug_s_quad = apn_ea_mugshot(s_quad)
                        same_mug_s_quad = db.query_functions({'mugshot': mug_s_quad})
                        same_mug_s_quad_id = set([t["id"] for t in same_mug_s_quad])
                        good_set = quad_id.intersection(same_mug_s_quad_id)
                        # We then test the EL-equivalence of ortho-derivatives
                        for j in good_set :

                            pi_s = ortho_derivative(s_quad)
                            pi_j = ortho_derivative(db[j]["sbox"].lut())
                            if are_linear_equivalent(pi_j,pi_s):

                                known_ccz_id.add(db[j]["ccz_id"])
                                known_id.update(set([fun["id"] for fun in db.query_functions({"ccz_id" : db[j]["ccz_id"]})]))
                                new_valve = False
                                break

                    # Otherwise, we just check EA equivalence  
                    else :
                        good_set = same_mug_s_id.intersection(non_ccz_quad_id)
                        for j in list(good_set):
                                if algebraic_degree(db[j]["sbox"].lut())==2 and algebraic_degree(s)==2:
                                    valve = are_ea_equivalent_from_vq(db[j]["sbox"].lut(),s)
                                else:
                                    valve = are_ea_equivalent(db[j]["sbox"].lut(),s)
                                if valve:
                                
                                    known_ccz_id.add(db[j]["ccz_id"])
                                    known_id.update(set([fun["id"] for fun in db.query_functions({"ccz_id" : db[j]["ccz_id"]})]))
                                    new_valve = False
                                    break
                if new_valve == True:
                    print("NEW FUNCTION")
                    print(s)
                    

    return(known_ccz_id,known_id)



# Helper to display the graph
def from_adj_list(G):
    n = len(G)
    g = nx.Graph()
    g.add_nodes_from(range(n))
    for u, nbrs in enumerate(G):
        for v in nbrs:
            g.add_edge(u, v)
    return g




n = 6
from time import *
with APNFunctions('apn6.db') as db:

    
    G = []
    all_apn = db.query_functions({})
    with Experiment("Computing the adjacency list of the CCZ-SN graph for n = 6"):

        ccz_dist_1, ea_dist_1 = set({}),set({})
        t= time()
        # For each known CCZ-class, we compute the CCZ-class it reaches with SN
        for i in range(14):
            section("Starting CCZ-class {}".format(i))
            
            ccz1 = db.query_functions({"ccz_id":i})
            pprint("The class contains {} functions".format(len(ccz1)))
            ccz_dist_1, ea_dist_1 = set({}),set({})
            # We get the SN for each EA-class in the CCZ-class and ass them to the known list
            for j in trange(len(ccz1)):
                f = ccz1[j]
                ccz_dist_1, ea_dist_1= at_distance_1_6(f,n,db,ccz_dist_1,ea_dist_1)
            
            pprint("The CCZ-class {} reached the following CCZ-classes:".format(i))
            print(list(ccz_dist_1))
            G.append(list(ccz_dist_1))

    
    # Plotting Graph. It's ugly ! At yout own risk
    """
    g = from_adj_list(G)
    pos = nx.spring_layout(g, k=1.2, iterations=50)  # k>default 0.1 spreads nodes
    nx.draw(g, pos, with_labels=True, node_size=50, font_size=5)
    plt.show()
    """
    section("Summary: the Adjacency List")
    for i in range(14):
        print("Class {} is connected to :".format(i),G[i])




###########
# Results #
###########

# This should be the adjacency list you obtain 
"""

at_dist_1 = [
    [0, 1],
    [0, 1],
    [2, 5, 6, 7, 8, 9, 10, 11, 12],
    [3, 6, 7, 9, 10, 11, 12],
    [4, 5, 8, 11],
    [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    [2, 3, 5, 6, 7, 8, 9, 10, 11, 12],
    [2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    [2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    [2, 3, 5, 6, 7, 8, 9, 10, 11, 12],
    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    [2, 3, 5, 6, 7, 8, 9, 10, 11, 12],
    [5, 7, 8, 9, 11, 13]
]

"""
