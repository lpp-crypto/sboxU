from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent
from matplotlib.widgets import Slider, RadioButtons

from sboxUv2.core import get_sbox
from sboxUv2.statistics import \
    ddt, differential_spectrum, ddt_coeff_probability, \
    lat, walsh_spectrum, absolute_walsh_spectrum, lat_coeff_probability_permutation, lat_coeff_probability_function, \
    bct, boomerang_spectrum, bct_coeff_probability, \
    fbct, fbct_spectrum \
    


# !SECTION! Interactive Tables



# !SUBSECTION! General function to generate interacting views

def table_interactive_view(
        table,
        title="Table",
        desc="v",
        cmap="viridis",
        vmin=0,
        vmax=40,
        with_sliders=True,
        with_cmap_choice=True,
):
    """Stops the execution flow and displays a new window containing a Pollock-style (in the sense of [C:BirPer15]) representation of a table.

    The table is assumed to be a 2-dimensional array containing integers. We associate a color to each value using the colormap `cmap`. Then, a 2D picture is generated representing the content of the table.

    A cross-hair cursor is displayed where the mouse, and the values of the coordinates are displayed under the picture. To describe the value at a given position, you need to specify the `desc` input. The cursor is handled using the TableCursor class.

    If `with_sliders` is set to True (the default case), the maximum value in the scale can be chosen with a slider.
    
    """
    fig, ax = plt.subplots()
    ax.set_title(title)
    im = ax.imshow(
        table,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    plt.subplots_adjust(right=0.85, left=0.25, bottom=0.2) # making space for all the widgets
    cbar = fig.colorbar(im,
                        cax=fig.add_axes([0.3, 0.02, 0.45, 0.02]),
                        orientation='horizontal')
    
    # setting up cursor
    cursor = TableCursor(ax,
                         lambda a,b : table[a][b] if a < len(table) and b < len(table[a]) else -1,
                         desc=desc)
    fig.canvas.mpl_connect('motion_notify_event', cursor.on_mouse_move)

    # Create sliders
    actual_max = 0
    actual_min = 0
    for row in table[1:]:
        actual_max = max(actual_max, max(row))
        actual_min = min(actual_max, min(row))
    if with_sliders:
        slider_max = Slider(
            ax=plt.axes([0.95, 0.2, 0.02, 0.5]),
            label='Scale max',
            valmin=actual_min * 0.8,
            valmax=actual_max * 1.2,
            valinit=actual_max,
            orientation="vertical",
        )
        slider_min = Slider(
            ax=plt.axes([0.90, 0.2, 0.02, 0.5]),
            label='Scale min',
            valmin=actual_min * 0.8,
            valmax=actual_max * 1.2,
            valinit=actual_min,
            orientation="vertical",
        )

        def update_from_sliders(val):
            im.set_clim(vmin=slider_min.val,
                        vmax=slider_max.val)
            fig.canvas.draw_idle()

        slider_max.on_changed(update_from_sliders)
        slider_min.on_changed(update_from_sliders)

    # Create cmap choice
    if with_cmap_choice:
        ax_radio = plt.axes([0.02, 0.3, 0.15, 0.4])
        radio = RadioButtons(
            ax_radio,
            labels=["viridis", "inferno", "cividis",
                    "coolwarm", "Spectral", "seismic", "berlin",
                    "Paired", "tab10", "tab20",
                    "gray"],
            active=0  # viridis is selected by default
        )

        def update_cmap(label):
            im.set_cmap(label)
            fig.canvas.draw_idle()

        radio.on_clicked(update_cmap)
    plt.show()
    


# !SUBSECTION! Specific instanciations 
    
def lat_interactive_view(
        s,
        cmap="coolwarm",
        vmin=None,
        vmax=40,
        absolute=True,
        show_only=None,
        with_sliders=True,
        with_cmap_choice=True
):
    sb = get_sbox(s)
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
    if show_only != None:
        t = [[t[a][b] if t[a][b] in show_only else vmax+2
              for b in range(0, len(t[a]))]
             for a in range(0, len(t))]
    table_interactive_view(t,
                           title="LAT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\sum_x (-1)^{ax +bS(x)}$",
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax,
                           with_sliders=with_sliders,
                           with_cmap_choice=with_cmap_choice
                           )
    
    
def ddt_interactive_view(
        s,
        cmap="coolwarm",
        vmin=0,
        vmax=10,
        show_only=None,
        with_sliders=True,
        with_cmap_choice=True
):
    sb = get_sbox(s)
    t = ddt(sb)
    if show_only != None:
        t = [[t[a][b] if t[a][b] in show_only else vmax+2
              for b in range(0, len(t[a]))]
             for a in range(0, len(t))]
    table_interactive_view(t,
                           title="DDT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\#\\{x, S(x+a)+S(x)=b\\}$",
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax,
                           with_sliders=with_sliders,
                           with_cmap_choice=with_cmap_choice
                           )
    
    
def bct_interactive_view(
        s,
        cmap="coolwarm",
        vmin=0,
        vmax=32,
        absolute=True,
        show_only=None,
        with_sliders=True,
        with_cmap_choice=True
):
    sb = get_sbox(s)
    t = bct(sb)
    if show_only != None:
        t = [[t[a][b] if t[a][b] in show_only else vmax+2
              for b in range(0, len(t[a]))]
             for a in range(0, len(t))]
    table_interactive_view(t,
                           title="BCT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\#\\{x, S^{-1}(S(x)+b) + S^{-1}(S(x+a)+b)=a\\}$",
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax,
                           with_sliders=with_sliders,
                           with_cmap_choice=with_cmap_choice
                           )

    
def fbct_interactive_view(
        s,
        cmap="coolwarm",
        vmin=0,
        vmax=16,
        absolute=True,
        show_only=None,
        with_sliders=True,
        with_cmap_choice=True
):
    sb = get_sbox(s)
    t = fbct(sb)
    if show_only != None:
        t = [[t[a][b] if t[a][b] in show_only else vmax+2
              for b in range(0, len(t[a]))]
             for a in range(0, len(t))]
    table_interactive_view(t,
                           title="F-BCT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\#\\{x, \\sum_{y \\in x+<a, b>}S(y)\\}$",
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax,
                           with_sliders=with_sliders,
                           with_cmap_choice=with_cmap_choice
                           )


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
        self.text = ax.text(0.2, -0.12, '', transform=ax.transAxes)
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

# !SUBSECTION! General function

def interactive_distribution_comparison(
        spec,
        in_length,
        out_length,
        expected_distrib,
        title="T",
        name="S",
        l_min=None,
        l_max=None,
        y_log_scale=True,
):
    # preprocessing arguments
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
    p.fill_between(abscissa, diff_min, diff_max, alpha=0.5, linewidth=0)
    # adding the metadata
    p.legend(shadow=True, fontsize=18)
    p.set_xlim([l_min, l_max])
    p.set_ylim([0.5, overall_max*1.2])
    p.yaxis.get_label().set_fontsize(18)
    p.xaxis.get_label().set_fontsize(18)
    p.tick_params(labelsize=12)
    p.grid(color="0.8")
    if y_log_scale:
        # if require_version(9):
        #     p.set_yscale("log", basey=2, nonposy="clip")
        # else:
        p.set_yscale("log", base=2, nonpositive="clip")
    plt.show()




# !SUBSECTION! Specific instanciations

def interactive_distribution_comparison_lat(s, y_log_scale=True):
    sb = get_sbox(s)
    if sb.is_invertible():
        interactive_distribution_comparison(
            absolute_walsh_spectrum(sb),
            sb.get_input_length(),
            sb.get_output_length(),
            lat_coeff_probability_permutation,
            title="LAT",
            name=sb.name().decode("UTF-8"),
            y_log_scale=y_log_scale
        )
    else:
        interactive_distribution_comparison(
            absolute_walsh_spectrum(sb),
            sb.get_input_length(),
            sb.get_output_length(),
            lat_coeff_probability_function,
            title="LAT",
            name=sb.name().decode("UTF-8"),
            y_log_scale=y_log_scale
        )
        
            

def interactive_distribution_comparison_ddt(s, y_log_scale=True):
    sb = get_sbox(s)
    interactive_distribution_comparison(
        differential_spectrum(sb),
        sb.get_input_length(),
        sb.get_output_length(),
        ddt_coeff_probability,
        title="DDT",
        name=sb.name().decode("UTF-8"),
        y_log_scale=y_log_scale
    )
            
 

def interactive_distribution_comparison_bct(s, y_log_scale=True):
    sb = get_sbox(s)
    interactive_distribution_comparison(
        boomerang_spectrum(sb),
        sb.get_input_length(),
        sb.get_output_length(),
        bct_coeff_probability,
        title="BCT",
        name=sb.name().decode("UTF-8"),
        y_log_scale=y_log_scale
    )
            
