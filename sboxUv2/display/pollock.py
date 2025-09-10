from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent

from sboxUv2.core import Sb
from sboxUv2.statistics import \
    ddt, differential_spectrum, ddt_coeff_probability, \
    lat, walsh_spectrum, absolute_walsh_spectrum, lat_coeff_probability_permutation, lat_coeff_probability_function, \
    bct, boomerang_spectrum, bct_coeff_probability \
    


# !SECTION! Interactive Tables


# !SUBSECTION! General function to generate interacting views

def table_interactive_view(
        table,
        title="Table",
        desc="v",
        cmap="coolwarm",
        vmin=0,
        vmax=20
):
    """Stops the execution flow and displays a new window containing a Pollock-style (in the sense of [C:BirPer15]) representation of a table.

    The table is assumed to be a 2-dimensional array containing integers. We associate a color to each value using the colormap `cmap`. Then, a 2D picture is generated representing the content of the table.

    A cross-hair cursor is displayed where the mouse, and the values of the coordinates are displayed under the picture. To describe the value at a given position, you need to specify the `desc` input. The cursor is handled using the TableCursor class.
    """
    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.imshow(
        table,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax
    )
    cursor = TableCursor(ax,
                         lambda a,b : table[a][b],
                         desc=desc)
    fig.canvas.mpl_connect('motion_notify_event', cursor.on_mouse_move)
    plt.show()


# !SUBSECTION! Specific instanciations 
    
def lat_interactive_view(s, cmap="coolwarm", vmin=None, vmax=40, absolute=True):
    sb = Sb(s)
    table = lat(sb)
    if absolute:
        t = [[abs(table[a][b]) for b in range(0, len(table[a]))]
             for a in range(0, len(table))]
        if vmin == None:
            vmin = 0
    else:
        t = table
        if vmin == None:
            vmin = -vmax
    table_interactive_view(t,
                           title="LAT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\sum_x (-1)^{ax +bS(x)}$",
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax)


    
    
def ddt_interactive_view(s, cmap="coolwarm", vmin=0, vmax=10, absolute=True):
    sb = Sb(s)
    table_interactive_view(ddt(sb),
                           title="DDT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\#\\{x, S(x+a)+S(x)=b\\}$",
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax)

    
    
def bct_interactive_view(s, cmap="coolwarm", vmin=0, vmax=32, absolute=True):
    sb = Sb(s)
    table_interactive_view(bct(sb),
                           title="BCT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\#\\{x, S^{-1}(S(x)+b) + S^{-1}(S(x+a)+b)=a\\}$",
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax)


# !SUBSECTION! The TableCursor class 
    
class TableCursor:
    """
    A cross hair cursor.

    Its implementation is largely taken from the matplotlib tutorial (see https://matplotlib.org/stable/gallery/event_handling/cursor_demo.html#sphx-glr-gallery-event-handling-cursor-demo-py).
    """
    def __init__(self, ax, values_function, desc):
        self.ax = ax
        self.horizontal_line = ax.axhline(color='k', lw=0.8, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.8, ls='--')
        # text location in axes coordinates
        self.text = ax.text(0.1, -0.1, '', transform=ax.transAxes)
        self.values_function = values_function
        self.desc = desc

        
    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    
    def on_mouse_move(self, event):
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            # update the line positions
            self.horizontal_line.set_ydata([y])
            self.vertical_line.set_xdata([x])
            a, b = int(round(y)), int(round(x))
            self.text.set_text('a=0x{:02x}   b=0x{:02x}   {}={:3d}'.format(
                a, b,
                self.desc,
                self.values_function(a, b)
            ))
            self.text.set_fontfamily("monospace")
            self.ax.figure.canvas.draw()



# !SECTION! Studying distributions of coefficients

def interactive_distribution_comparison(
        spec,
        in_length,
        out_length,
        expected_distrib,
        title="T",
        name="S",
        l_min=None,
        l_max=None,
        x_log_scale=False,
        y_log_scale=True
):
    spectra = {}
    spectra[name] = defaultdict(float)
    absolute_spec = defaultdict(int)
    for k in spec.keys():
        spectra[name][abs(k)] += spec[k]
    if l_min == None:
        l_min = 0
    # computing expected probabilities
    spectra["Expected"] = defaultdict(float)
    finished = False
    c = l_min
    while not finished:
        p = expected_distrib(in_length, out_length, c)
        expected_card = p*2**out_length*(2**in_length-1)
        if expected_card > 2**-4:
            spectra["Expected"][c] = expected_card
        elif p != 0:
            finished = True
        c += 2
    l_max = max(c, max(spectra[name].keys()) + 2)

    # plotting data itself
    fig, p = plt.subplots()
    p.set_xlabel("c")
    p.set_ylabel("$\\# \\{(a, b), |" + title + "[a,b]| = c\\}$")
    overall_max = 9
    for w in spectra.keys():
        abscissa = []
        ordenna = []
        for c in sorted(spectra[w].keys()):
            if c >= l_min and c <= l_max:
                abscissa.append(c)
                v = float(spectra[w][c])
                ordenna.append(v)
                overall_max = max(overall_max, v)
        p.plot(abscissa,
                   ordenna,
                   marker="x" if w == "Expected" else "o",
                   linestyle="-" if w == "Expected" else "",
                   markersize=4,
                   label=w)
    # plotting the difference
    abscissa = []
    diff_min = []
    diff_max = []
    for c in range(l_min, l_max+1):
        if spectra[name][c] != 0 or spectra["Expected"][c] != 0:
            abscissa.append(float(c))
            diff_min.append(float(min(spectra[name][c], spectra["Expected"][c])))
            diff_max.append(float(max(spectra[name][c], spectra["Expected"][c])))
    print(abscissa)
    print(diff_min)
    print(diff_max)
    p.fill_between(abscissa, diff_min, diff_max, alpha=0.5, linewidth=0)
    # adding the metadata
    p.legend(shadow=True, fontsize=18)
    p.set_xlim([l_min, l_max])
    p.set_ylim([0.5, overall_max*1.2])
    p.yaxis.get_label().set_fontsize(18)
    p.xaxis.get_label().set_fontsize(18)
    p.tick_params(labelsize=12)
    p.grid(color="0.8")
    if x_log_scale:
        # if require_version(9):
            # p.set_xscale("log", basex=2, nonposx="clip")
        # else:
        p.set_xscale("log", base=2, nonpositive="clip")
    if y_log_scale:
        # if require_version(9):
        #     p.set_yscale("log", basey=2, nonposy="clip")
        # else:
        p.set_yscale("log", base=2, nonpositive="clip")
    plt.show()


def interactive_distribution_comparison_lat(s):
    sb = Sb(s)
    if sb.is_invertible():
        interactive_distribution_comparison(
            absolute_walsh_spectrum(sb),
            sb.get_input_length(),
            sb.get_output_length(),
            lat_coeff_probability_permutation,
            title="LAT",
            name=sb.name().decode("UTF-8")
        )
    else:
        interactive_distribution_comparison(
            absolute_walsh_spectrum(sb),
            sb.get_input_length(),
            sb.get_output_length(),
            lat_coeff_probability_function,
            title="LAT",
            name=sb.name().decode("UTF-8")
        )
        
            
