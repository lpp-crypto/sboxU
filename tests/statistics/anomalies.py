from sage.all import *
from sboxUv2 import *
from collections import defaultdict

# ----

n = 8
n_tested = 2**13
with Experiment("Testing distribution of maximum values"):

# ----

    section("LAT")

# ----

    subsection("Experimental")
    dis = defaultdict(int)
    for t in range(0, n_tested):
        dis[linearity(random_permutation_S_box(n))] += 1
    for k in sorted(dis.keys()):
        if dis[k] != 0:
            print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))

# ----

    subsection("Theoretical")
    dis = expected_linearity_distribution_permutation(n, n)
    for k in sorted(dis.keys()):
        if dis[k] != 0.0:
            print("{:3d}: {:.4f}".format(k, dis[k]))

# ----

    section("DDT")

# ----

    subsection("Experimental")
    dis = defaultdict(int)
    for t in range(0, n_tested):
        dis[differential_uniformity(random_permutation_S_box(n))] += 1
    for k in sorted(dis.keys()):
        if dis[k] != 0:
            print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
    
    subsection("Theoretical")
    dis = expected_differential_uniformity_distribution_permutation(n, n)
    for k in sorted(dis.keys()):
        if dis[k] != 0.0:
            print("{:3d}: {:.4f}".format(k, dis[k]))

# ----

    section("BCT")

# ----

    subsection("Experimental")
    dis = defaultdict(int)
    for t in range(0, n_tested):
        dis[boomerang_uniformity(random_permutation_S_box(n))] += 1
    for k in sorted(dis.keys()):
        if dis[k] != 0:
            print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
    
    subsection("Theoretical")
    dis = expected_boomerang_uniformity_distribution_permutation(n, n)
    for k in sorted(dis.keys()):
        if dis[k] != 0.0:
            print("{:3d}: {:.4f}".format(k, dis[k]))

# ----

