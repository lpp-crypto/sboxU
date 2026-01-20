from sage.all import *
from sboxUv2 import *
import time
import numpy as np
import matplotlib.pyplot as plt

###########
# Toolbox #
###########

def time_f(n, m, n_repeat):
    """Return median runtime (seconds) of f(n, m) over `repeats` runs."""
    times = []
    for _ in range(n_repeat):
        t0 = time.perf_counter()
        non_trivial_sn(0,Sb(X**3),n, m)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return np.mean(times)

def find_best_nm(n,s,n_repeat,
    first_batch_start,
    first_batch_end,
    first_batch_step,
    second_batch_start, 
    second_batch_end,
    second_batch_step,
    show_plot=True,
):
    """
    Search over n_range x m_range, return (best_n, best_m, time_matrix).
    n_range, m_range: sequences (e.g. range(...) or list of ints)
    """
    n_range = range(first_batch_start,first_batch_end,first_batch_step)
    m_range = range(second_batch_start,second_batch_end,second_batch_step)

    with Experiment("Looking for the best parameters for SN computation"):

        section("Parameters")
        print("Testing over {} bits with APN function :".format(n))
        print(s.lut())
        print()
        print("For the first batch, testing from {} to {} with step {}".format(first_batch_start,
    first_batch_end,
    first_batch_step))
        print("For the first batch, testing from {} to {} with step {}".format(second_batch_start, 
        second_batch_end,
        second_batch_step))

        n_vals = list(n_range)
        m_vals = list(m_range)
        T = np.zeros((len(n_vals), len(m_vals)), dtype=float)
        section("Starting Tests")
        print("It might take some time, start with a greater step and smaller range for your first tests")
        pprint("[bold]Number of Evaluations: {}[/bold]".format(n_repeat*len(n_vals)*len(m_vals)))
        print()
        for i, n in enumerate(n_vals):
            for j, m in enumerate(m_vals):
                T[i, j] = time_f(n, m, n_repeat)

        # find minimum time
        idx = np.unravel_index(np.argmin(T), T.shape)
        best_n = n_vals[idx[0]]
        best_m = m_vals[idx[1]]
        best_time = T[idx]
        pprint("[bold]Best Value for First Batch: {}[/bold]".format(best_n))
        pprint("[bold]Best Value for Second Batch: {}[/bold]".format(best_m))
        pprint("[green]Best Time: {}[/green]".format(best_time))
        print()
        print("Now Showing the Time Distribution")
        pprint("[red]CLOSE THE GRAPH TO EXIT[/red]")
        if show_plot:
            N, M = np.meshgrid(m_vals, n_vals)  # X: m, Y: n
            plt.figure(figsize=(6, 4))
            cp = plt.contourf(N, M, T, levels=20, cmap="viridis")
            plt.colorbar(cp, label="Runtime (s)")
            plt.scatter([best_m], [best_n], color="red", label="Best (n,m)")
            plt.xlabel("m")
            plt.ylabel("n")
            plt.title("Runtime of f(n, m)")
            plt.legend()
            plt.tight_layout()
            plt.show()

    return best_n, best_m, best_time

##############
#### MAIN ####
##############

if __name__ == "__main__":

    ##############
    # Parameters #
    ##############

    # Choose n
    n = 8
    # Choose an APN function
    field = GF(2**n)
    X = PolynomialRing(field, "X").gen()
    s = Sb(X**3)
    # Choose a range of size for the first batch 
    first_batch_start = 625
    first_batch_end = 640
    first_batch_step = 5
    # Choose a range of size for the first batch 
    second_batch_start = 50
    second_batch_end = 90
    second_batch_step = 5
    # Choose a number of time to repeat for each parameters
    n_repeat = 100

    best_n, best_m, T = find_best_nm(n,s,n_repeat,
                                    first_batch_start,
                                    first_batch_end,
                                    first_batch_step,
                                    second_batch_start, 
                                    second_batch_end,
                                    second_batch_step)
