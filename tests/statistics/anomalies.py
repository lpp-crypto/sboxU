from sage.all import *
from sboxUv2 import *
from collections import defaultdict

# ----

from sage.crypto.sboxes import sboxes

# ----

with Experiment("Checking some anomalies"):

# ----

    section("Skipjack")
    skipjack = Sb(sboxes["Skipjack"])
    for table in ["LAT", "BCT", "DDT"]:
        print("{} {:10.5f} {:10.5f}".format(
            table,
            table_anomaly(skipjack, table),
            table_negative_anomaly(skipjack, table)
        ))

# ----

    interactive_distribution_comparison_lat(skipjack)

# ----

    section("Other S-boxes")
    for k in ["AES", "Fantomas", "Kuznyechik"]:
        s = Sb(sboxes[k])
        subsection(k)
        for table in ["LAT", "BCT", "DDT"]:
            print("{} {:10.5f} {:10.5f}".format(
                table,
                table_anomaly(s, table),
                table_negative_anomaly(s, table)
            ))

# ----

with Experiment("Testing distribution of maximum values"):

# ----

    n = 8
    n_tested = 2**13

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

