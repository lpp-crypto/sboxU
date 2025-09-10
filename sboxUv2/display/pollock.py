import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent

from sboxUv2.core import Sb
from sboxUv2.statistics import ddt, lat, bct


def table_interactive_view(
        table,
        title="Table",
        desc="v",
        cmap="coolwarm",
        vmin=0,
        vmax=20
):
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

    
    
def lat_interactive_view(s, cmap="coolwarm", vmin=0, vmax=40, absolute=True):
    sb = Sb(s)
    table = lat(sb)
    if absolute:
        t = [[abs(table[a][b]) for b in range(0, len(table[a]))]
             for a in range(0, len(table))]
    else:
        t = table
    table_interactive_view(t,
                           title="LAT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\sum_x (-1)^{ax +bS(x)}$",
                           vmin=vmin,
                           vmax=vmax)


    
    
def ddt_interactive_view(s, cmap="coolwarm", vmin=0, vmax=10, absolute=True):
    sb = Sb(s)
    table_interactive_view(ddt(sb),
                           title="DDT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\#\\{x, S(x+a)+S(x)=b\\}$",
                           vmin=vmin,
                           vmax=vmax)

    
    
def bct_interactive_view(s, cmap="coolwarm", vmin=0, vmax=32, absolute=True):
    sb = Sb(s)
    table_interactive_view(bct(sb),
                           title="BCT of {}".format(sb.name().decode("UTF-8")),
                           desc="$\\#\\{x, S^{-1}(S(x)+b) + S^{-1}(S(x+a)+b)=a\\}$",
                           vmin=vmin,
                           vmax=vmax)


    
class TableCursor:
    """
    A cross hair cursor.
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


