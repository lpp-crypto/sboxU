# ~/~ begin <<anomalies.md#maximums.py>>[init]
from sage.all import *
from sboxU import *
from collections import defaultdict
# ~/~ end
# ~/~ begin <<anomalies.md#maximums.py>>[1]
n = 8
n_tested = 2**13
# ~/~ end
# ~/~ begin <<anomalies.md#maximums.py>>[2]
dis = defaultdict(int)
print("---- experimental")
for t in range(0, n_tested):
    dis[linearity(random_permutation_S_box(n))] += 1
for k in sorted(dis.keys()):
    if dis[k] != 0:
        print("{:3d}: {:.4f}".format(k, float(dis[k])/n_tested))
print("---- theoretical")
dis = expected_linearity_distribution_permutation(n, n)
for k in sorted(dis.keys()):
    if dis[k] != 0.0:
        print("{:3d}: {:.4f}".format(k, dis[k]))
# ~/~ end
