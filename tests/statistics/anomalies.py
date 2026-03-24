import sys
from sage.all import *
from sboxU import *
# --- { 
from sage.crypto.sboxes import sboxes
from collections import defaultdict
# --- } 
def main_test():
    with Experiment(' Computing some anomalies'):
        section(' Computing some anomalies')
        subsection(' Skipjack')
        # --- { 
        skipjack = get_sbox(sboxes["Skipjack"])
        for table in ["LAT", "BCT", "DDT"]:
            print("{} {:10.5f} {:10.5f}".format(
                table,
                table_anomaly(skipjack, table),
                table_negative_anomaly(skipjack, table)
            ))
        # --- } 
        # --- { 
        interactive_distribution_comparison_lat(skipjack)
        # --- } 
        subsection(' Other S-boxes')
        # --- { 
        for k in ["AES", "Fantomas", "Kuznyechik"]:
            s = get_sbox(sboxes[k])
            print(k)
            for table in ["LAT", "BCT", "DDT"]:
                print("{} {:10.5f} {:10.5f}".format(
                    table,
                    table_anomaly(s, table),
                    table_negative_anomaly(s, table)
                ))
        # --- } 
        section(' Looking at the maximum value')
        # --- { 
        n = 8
        n_tested = 2**11
        # --- } 
        subsection(' LAT ')
        print(' Experimental ')
        # --- { 
        dis = defaultdict(int)
        for t in range(0, n_tested):
            dis[linearity(random_permutation_S_box(n))] += 1
        for k in sorted(dis.keys()):
            if dis[k] != 0:
                print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
        # --- } 
        print(' Theoretical')
        # --- { 
        dis = expected_linearity_distribution_permutation(n, n)
        for k in sorted(dis.keys()):
            if dis[k] != 0.0:
                print("{:3d}: {:.4f}".format(k, dis[k]))
        # --- } 
        subsection(' DDT')
        # --- { 
        print("Experimental")
        dis = defaultdict(int)
        for t in range(0, n_tested):
            dis[differential_uniformity(random_permutation_S_box(n))] += 1
        for k in sorted(dis.keys()):
            if dis[k] != 0:
                print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
        print("Theoretical")
        dis = expected_differential_uniformity_distribution_permutation(n, n)
        for k in sorted(dis.keys()):
            if dis[k] != 0.0:
                print("{:3d}: {:.4f}".format(k, dis[k]))
        # --- } 
        subsection(' BCT')
        # --- { 
        print("Experimental")
        dis = defaultdict(int)
        for t in range(0, n_tested):
            dis[boomerang_uniformity(random_permutation_S_box(n))] += 1
        for k in sorted(dis.keys()):
            if dis[k] != 0:
                print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
        print("Theoretical")
        dis = expected_boomerang_uniformity_distribution_permutation(n, n)
        for k in sorted(dis.keys()):
            if dis[k] != 0.0:
                print("{:3d}: {:.4f}".format(k, dis[k]))
        fail("bwooo")
        # --- } 
        subsection(' Comments')
    return exit_code()


if __name__ == '__main__':    sys.exit(main_test())
