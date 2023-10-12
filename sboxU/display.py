#!/usr/bin/sage
# Time-stamp: <2023-10-12 17:09:16 lperrin>

import matplotlib
import matplotlib.pyplot as plt
from .diff_lin import *
from math import log
from collections import defaultdict

# source for color sequence: https://matplotlib.org/users/dflt_style_changes.html

COLOR_SEQUENCE = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]
DEFAULT_FIG_X, DEFAULT_FIG_Y = 15, 15


# !SECTION! Display in the console 

def pretty_spectrum(d, absolute=False):
    """Returns a line containing a pretty representation of the
    dictionnary d.

    """
    if len(d.keys()) == 0:
        return "{}"

    if absolute == False:
        printed_dict = d
    else:
        printed_dict = defaultdict(int)
        for k in d.keys():
            printed_dict[abs(k)] += d[k]
    line = "{"
    for k in sorted(printed_dict.keys()):
        line += "{}: {}, ".format(k, printed_dict[k])
    return line[:-2] + "}"


def pretty_vector(v, template="{:2x},"):
    """Returns a string containing the representation of the integers in v
    using the template given (defaults to a simple decimal
    representation).

    """
    if len(v) == 0:
        return "[]"
    line = "["
    for x in v:
        line += template.format(x)
    return line[:-1] + "]"

    
def pretty_lagrange(s, G):
    poly_ring = PolynomialRing(G, "y")
    p = poly_ring.lagrange_polynomial([(G.fetch_int(i), G.fetch_int(s[i])) for i in range(0, len(s))])
    result = ""
    for i, k in enumerate(p):
        if k == 1:
            result += "X^{:d} + ".format(i)
        elif k != 0:
            result += "{:x}*X^{:d} + ".format(k.integer_representation(), i)
    return result[:-2]


# !SECTION! Graph generation

def plot_table_averages(l,
                        file_name="avg",
                        rows=True,
                        cols=True,
                        col_rows=None,
                        col_cols=None):
    if rows and cols:
        col_rows, col_cols = COLOR_SEQUENCE[0], COLOR_SEQUENCE[1]
    elif rows:
        col_rows = COLOR_SEQUENCE[0]
    elif cols:
        col_cols = COLOR_SEQUENCE[0]
    else:
        raise "At least rows or cols should be set to True!"
    fig, p = plt.subplots(figsize=(DEFAULT_FIG_X, DEFAULT_FIG_Y))
    p.set_xlabel('row/column index')
    p.set_ylabel('Average of the absolute value')
    # rows
    if rows:
        avgs_a = []
        for a in range(1, len(l)):
            avg = 0.0
            for w in l[a]:
                avg += abs(float(w))
            avg = avg/len(l[a])
            avgs_a.append(avg)
        p.plot(range(1, len(avgs_a)+1),
               avgs_a,
               marker="o",
               color=col_rows,
               linestyle="-",
               markersize=2,
               label="Rows")
    # columns
    if cols:
        avgs_b = []
        for b in range(1, len(l)):
            avg = 0.0
            col = [l[a][b] for a in range(0, len(l))]
            for w in col:
                avg += abs(float(w))
            avg = avg/len(col)
            avgs_b.append(avg)
        p.plot(range(1, len(avgs_b)+1),
               avgs_b,
               color=col_cols,
               marker="^",
               linestyle="-",
               markersize=2,
               label="Columns")
    # finalizing graph
    legend = p.legend(loc='upper right', shadow=True)
    p.set_xlim([0, len(l)])
    fig.savefig("{}.png".format(file_name))
    plt.close()


def plot_table_variances(l,
                        file_name="var",
                        rows=True,
                        cols=True,
                        col_rows=None,
                        col_cols=None):
    if rows and cols:
        col_rows, col_cols = COLOR_SEQUENCE[0], COLOR_SEQUENCE[1]
    elif rows:
        col_rows = COLOR_SEQUENCE[0]
    elif cols:
        col_cols = COLOR_SEQUENCE[0]
    else:
        raise "At least rows or cols should be set to True!"
    fig, p = plt.subplots(figsize=(DEFAULT_FIG_X, DEFAULT_FIG_Y))
    p.set_xlabel('row/column index')
    p.set_ylabel('Variance of the absolute value')
    # rows
    if rows:
        variances_a = []
        for a in range(1, len(l)):
            avg = 0.0
            for w in l[a]:
                avg += abs(float(w))
            avg = avg/len(l[a])
            v = 0.0
            for w in l[a]:
                v += (abs(w) - avg)**2
            v = v / float(len(l[a]))
            variances_a.append(v)
        p.plot(range(1, len(variances_a)+1),
               variances_a,
               marker="o",
               color=col_rows,
               linestyle="-",
               markersize=2,
               label="Rows")
    # columns
    if cols:
        variances_b = []
        for b in range(1, len(l)):
            avg = 0.0
            col = [l[a][b] for a in range(0, len(l))]
            for w in col:
                avg += abs(float(w))
            avg = avg/len(col)
            v = 0.0
            for w in col:
                v += (abs(w) - avg)**2
            v = v / float(len(col))
            variances_b.append(v)
        p.plot(range(1, len(variances_b)+1),
               variances_b,
               color=col_cols,
               marker="o",
               linestyle="-",
               markersize=2,
               label="Columns")
    # finalizing graph
    p.legend(loc='upper right', shadow=True, fontsize=18)
    p.yaxis.get_label().set_fontsize(18)
    p.xaxis.get_label().set_fontsize(18)
    p.tick_params(labelsize=12)
    p.set_xlim([0, len(l)])
    n = int(round(log(len(l), 2)))
    if n > 4:
        tick_space = int(2**(n-4))
    else:
        tick_space = 1
    plt.xticks(list(range(0, 2**n, tick_space)))
    p.grid(color="0.8")
    fig.savefig("{}.png".format(file_name))
    plt.close()


def plot_differential(dict_s,
                      file_name="differential",
                      with_random_permutation=True,
                      with_random_function=False,
                      u_max=12,
                      x_log_scale=False,
                      y_log_scale=True):
    # distribution is the same for random functions and permutations
    with_random = with_random_function or with_random_permutation
    # coeff probabilities for s
    spectra = {}
    for func_name in dict_s.keys():
        spectra[func_name] = defaultdict(float)
        s = dict_s[func_name]
        diff_spec = differential_spectrum(s)
        for k in diff_spec.keys():
            spectra[func_name][k] = float(diff_spec[k]) / (len(s) * (len(s) - 1))
    # coeff probabilities for random function
    n = float(log(len(s), 2))
    if with_random:
        spectra["Random Permutation"] = {
            c: ddt_coeff_probability(n, n, c)
            for c in range(0, u_max+1, 2)
        }
    # plotting
    abscissa = range(0, u_max+1, 2)
    fig, p = plt.subplots(figsize=(DEFAULT_FIG_X, DEFAULT_FIG_Y))
    p.set_xlabel('DDT coefficients')
    p.set_ylabel('Number of occurrences')
    if x_log_scale:
        p.set_xscale("log", base=2, nonposx="clip")
    if y_log_scale:
        p.set_yscale("log", base=2, nonposy="clip")
    color_index = 0
    for w in spectra.keys():
        ordenna = []
        for c in range(0, u_max+1, 2):
            if c in spectra[w].keys():
                ordenna.append(float(spectra[w][c]))
            else:
                ordenna.append(0.0)
        p.plot(abscissa,
               ordenna,
               color=COLOR_SEQUENCE[color_index],
               marker="o",
               linestyle="-",
               markersize=2,
               label=w)
        color_index +=1
        legend = p.legend(loc='upper right', shadow=True)
    p.set_xlim([0, u_max])
    fig.savefig("{}.png".format(file_name))
    plt.close()
    
        

def plot_linear(dict_s,
                file_name="linear",
                with_random_permutation=True,
                with_random_function=False,
                l_min=0,
                l_max=64,
                x_log_scale=False,
                y_log_scale=True):
    # coeff probabilities for s
    spectra = {}
    for func_name in dict_s.keys():
        spectra[func_name] = defaultdict(float)
        s = dict_s[func_name]
        walsh_spec = walsh_spectrum(s)
        for k in walsh_spec.keys():
            if k == 0:
                spectra[func_name][k] = float(walsh_spec[k]) / (len(s) * (len(s) - 1))
            else:
                spectra[func_name][abs(k)] += float(walsh_spec[k]) / (len(s) * (len(s) - 1))
    # coeff probabilities for random functions
    n = float(log(len(s), 2))
    if with_random_function:
        spectra["Random Function"] = {
            c: lat_coeff_probability_function(n, n, c)
            for c in range(l_min, l_max+1, 4)
        }
    n = float(log(len(s), 2))
    if with_random_permutation:
        spectra["Random Permutation"] = {
            c: lat_coeff_probability_permutation(n, n, c)
            for c in range(l_min, l_max+1, 4)
        }
    # plotting
    abscissa = range(l_min, l_max+1, 4)
    fig, p = plt.subplots(figsize=(DEFAULT_FIG_X,DEFAULT_FIG_Y))
    p.set_xlabel('abs(LAT coefficients)')
    p.set_ylabel('Number of occurrences')
    color_index = 0
    for w in spectra.keys():
        ordenna = []
        for c in range(l_min, l_max+1, 4):
            if c in spectra[w].keys():
                ordenna.append(float(spectra[w][c]))
            else:
                ordenna.append(0.0)
        p.plot(abscissa,
               ordenna,
               color=COLOR_SEQUENCE[color_index],
               marker="o",
               linestyle="-",
               markersize=2,
               label=w)
        color_index +=1
        legend = p.legend(loc='upper right', shadow=True)
    p.set_xlim([l_min, l_max])
    if x_log_scale:
        p.set_xscale("log", base=2, nonposx="clip")
    if y_log_scale:
        p.set_yscale("log", base=2, nonposy="clip")
    fig.savefig("{}.png".format(file_name))
    plt.close()

        

def plot_statistical(spec,
                     n=None,
                     file_name="statistical",
                     expected_distrib=None,
                     l_min=None,
                     l_max=None,
                     x_log_scale=False,
                     y_log_scale=True):
    spectra = {}
    spectra["F"] = defaultdict(float)
    absolute_spec = defaultdict(int)
    for k in spec.keys():
        spectra["F"][abs(k)] += spec[k]
    if l_min == None:
        l_min = 0
    if expected_distrib != None:
        if n == None:
            raise Exception("in plot_statistical: `n` must be specified!")
        spectra["Expected"] = {}
        finished = False
        c = l_min
        while not finished:
            p = expected_distrib(n, n, c)
            expected_card = p*2**n*(2**n-1)
            if expected_card > 2**-4:
                spectra["Expected"][c] = expected_card
            elif p != 0:
                finished = True
            c += 2
        l_max = max(c, max(spectra["F"].keys()) + 2)
    else:
        l_max = max(spectra["F"].keys())+2

    # plotting
    fig, p = plt.subplots(figsize=(DEFAULT_FIG_X, DEFAULT_FIG_Y))
    p.set_xlabel("c")
    p.set_ylabel("# {(a, b), |T[a,b]| = c}")
    color_index = 0
    for w in spectra.keys():
        abscissa = []
        ordenna = []
        for c in sorted(spectra[w].keys()):
            if c >= l_min and c <= l_max:
                abscissa.append(c)
                ordenna.append(float(spectra[w][c]))
        if w == "Expected":
            p.plot(abscissa,
                   ordenna,
                   color=COLOR_SEQUENCE[color_index],
                   marker="x",
                   linestyle="-",
                   markersize=4,
                   label=w)
        else:
            p.plot(abscissa,
                   ordenna,
                   color=COLOR_SEQUENCE[color_index],
                   marker="o",
                   linestyle="",
                   markersize=8,
                   label=w)
        color_index +=1
    p.legend(shadow=True, fontsize=18)
    p.set_xlim([l_min, l_max])
    p.yaxis.get_label().set_fontsize(18)
    p.xaxis.get_label().set_fontsize(18)
    p.tick_params(labelsize=12)
    p.grid(color="0.8")
    if x_log_scale:
        p.set_xscale("log", base=2, nonpositive="clip")
    if y_log_scale:
        p.set_yscale("log", base=2, nonpositive="clip")
    fig.savefig("{}.png".format(file_name))
    plt.close()



# !SUBSECTION! Properties of rows and columns

def plot_statistical_by_rows(t, 
                             n=None,
                             file_name="statistical",
                             expected_distrib=None,
                             l_min=None,
                             l_max=None,
                             x_log_scale=False,
                             y_log_scale=True):
    spectra = {}
    if n == None:
        n = int(round(log(len(t), 2)))
    min_c = 2**n
    max_c = 0
    for a in range(1, 2**n):
        spectra[a] = defaultdict(int)
        for b in range(1, 2**n):
            c = abs(t[a][b])
            spectra[a][c] += 1
            max_c = max(c, max_c)
            min_c = min(c, min_c)
    if l_min == None:
        l_min = 0
    if expected_distrib != None:
        if n == None:
            raise Exception("in plot_statistical: `n` must be specified!")
        spectra[-1] = {}
        finished = False
        c = l_min
        while not finished:
            p = expected_distrib(n, n, c)
            expected_card = p*2**n
            if expected_card > 0.8:
                spectra[-1][c] = expected_card
            elif p > 0:
                finished = True
            c += 2
        if l_max == None:
            l_max = max(c, max_c+2)

    # plotting
    fig, p = plt.subplots(figsize=(DEFAULT_FIG_X, DEFAULT_FIG_Y))
    p.set_xlabel("c")
    p.set_ylabel("# {b, |T[a,b]| = c}")
    color_index = 0
    local_color_sequence = matplotlib.cm.rainbow(range(0, 2**n+1))
    for w in reversed(sorted(spectra.keys())):
        abscissa = []
        ordenna = []
        for c in sorted(spectra[w].keys()):
            if c >= l_min and c <= l_max:
                abscissa.append(c)
                ordenna.append(float(spectra[w][c]))
        if w == -1:
            p.plot(abscissa,
                   ordenna,
                   color="black",
                   marker="x",
                   linestyle="-",
                   markersize=4,
                   label=w)
        else:
            p.plot(abscissa,
                   ordenna,
                   color=local_color_sequence[color_index],
                   marker="o",
                   linestyle="-",
                   markersize=2,
                   alpha=0.2,
                   label=w)
        color_index +=1
    p.set_xlim([l_min, l_max])
    p.yaxis.get_label().set_fontsize(18)
    p.xaxis.get_label().set_fontsize(18)
    p.tick_params(labelsize=12)
    p.grid(color="0.8")
    if x_log_scale:
        p.set_xscale("log", base=2, nonpositive="clip")
    if y_log_scale:
        p.set_yscale("log", base=2, nonpositive="clip")
    fig.savefig("{}.png".format(file_name))
    plt.close()

    


# !SECTION! Jackson Pollock


def xor_texture(mat):
    """Returns a list of lists `t` such that

    t[i][j] = # { (x,y) | mat[x+i][y+j] = mat[x][y] } ,

    where "+" denotes the bitwise XOR.

    """
    coeff_coords = defaultdict(list)
    for a in range(0, len(mat)):
        for b in range(0, len(mat[a])):
            coeff_coords[abs(mat[a][b])].append([a, b])
    result = [[0 for x in range(0, len(mat[0]))]
              for y in range(0, len(mat))]
    for c in coeff_coords.keys():
        for coord0, coord1 in itertools.combinations(coeff_coords[c], r=2):
            i = oplus(coord0[0], coord1[0])
            j = oplus(coord0[1], coord1[1])
            result[i][j] += 2
    return result


def add_texture(mat):
    """Returns a list of lists `t` such that

    t[i][j] = # { (x,y) | mat[x+i][y+j] = mat[x][y] } ,

    where "+" denotes the addition modulo the number of lines/rows.

    """
    coeff_coords = defaultdict(list)
    for a in range(0, len(mat)):
        for b in range(0, len(mat[a])):
            coeff_coords[abs(mat[a][b])].append([a, b])
    result = [[0 for x in range(0, len(mat[0]))]
              for y in range(0, len(mat))]
    for c in coeff_coords.keys():
        for coord0, coord1 in itertools.product(coeff_coords[c], repeat=2):
            if coord0 != coord1:
                i = (coord0[0] - coord1[0]) % len(mat)
                j = (coord0[1] - coord1[1]) % len(mat)
                result[i][j] += 1
                result[(- i) % len(mat)][(-j) % len(mat)] += 1
    return result

    
def save_pollock(mat,
                 color_scheme="CMRmap_r",
                 name="pollock",
                 vmin=0,
                 vmax=20,
                 folder=None,
                 frame=True,
                 visible_axes=False,
                 colorbar=False,
                 file_type="png",
                 modifier_func=abs,
                 figsize=DEFAULT_FIG_X,
                 title=None):
    fig, p = plt.subplots(figsize=(figsize,figsize))
    abs_mat = [[modifier_func(mat[i][j]) for j in range(0, len(mat[0]))]
               for i in range(0, len(mat))]
    axes = p.imshow(
        abs_mat,
        interpolation="None",
        cmap=plt.cm.get_cmap(color_scheme, 100),
        vmin=vmin,
        vmax=vmax
    )
    p.set_aspect('equal')
    p.get_xaxis().set_visible(visible_axes)
    p.get_yaxis().set_visible(visible_axes)
    p.patch.set_alpha(0)
    p.set_frame_on(frame)
    if title != None:
        plt.title(title, fontsize=20)
    if colorbar:
        cax = fig.add_axes([p.get_position().x1+0.01, 
                            p.get_position().y0,
                            0.02, 
                            p.get_position().height])
        plt.colorbar(axes, cax=cax)
    if folder == None:
        name_base = "{}."+file_type
    else:
        name_base = folder + "/{}." + file_type
    fig.savefig(name_base.format(name))
    plt.close()

