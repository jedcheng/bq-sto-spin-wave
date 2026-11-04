"""Overlap checker for the figure scripts:  python check_overlap.py figX.py
Runs the script, and just before each savefig lists
  * text artists (labels, tick labels, titles, annotations, legend entries) whose boxes overlap each other,
  * text that sits on the plotting area of a different axes,
Tick labels of ticks outside the axis range (kept by matplotlib but not drawn), hidden text and empty strings are
ignored.  Text past the figure edge is not reported: the figures are saved with bbox_inches="tight".  The figure is still saved as usual."""
import sys, runpy, itertools
import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.text import Text

PAD = 1.0          # px shrink of every box, so touching boxes are not reported


def texts(fig):
    out, ticklab, drawn = [], set(), set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tk in axis.get_major_ticks() + axis.get_minor_ticks(): ticklab |= {id(tk.label1), id(tk.label2)}
            if not axis.get_visible() or not ax.axison: continue              # twinx x axis, axis("off") images
            for tk in axis._update_ticks():
                if tk.label1.get_visible(): drawn.add(id(tk.label1))
                if tk.label2.get_visible(): drawn.add(id(tk.label2))
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip(): continue
        if id(t) in ticklab and id(t) not in drawn: continue
        ax = t.axes
        if ax is not None and not ax.get_visible(): continue
        try: bb = t.get_window_extent()
        except Exception: continue
        if bb.width < 1 or bb.height < 1: continue
        out.append((t, bb.expanded(1, 1).padded(-PAD), ax))
    return out


def name(t):
    s = t.get_text().replace("\n", " | ")
    return f"'{s[:40]}'"


def check(fig):
    fig.canvas.draw()
    T = texts(fig); n = 0
    for (a, ba, axa), (b, bb, axb) in itertools.combinations(T, 2):
        if ba.overlaps(bb):
            print(f"  text/text   {name(a)}  x  {name(b)}"); n += 1
    axes = [ax for ax in fig.axes if ax.get_visible() and ax.axison]
    for t, bt, at in T:
        for ax in axes:
            if ax is at: continue
            if at is not None and at._twinned_axes.joined(at, ax): continue      # twinx partner
            if bt.overlaps(ax.bbox.padded(-PAD)) and ax.get_label() != "<colorbar>":
                print(f"  text/axes   {name(t)}  on axes {ax.get_label() or id(ax)} ({ax.get_xlabel() or ax.get_ylabel() or ax.get_title()})"); n += 1
    print(f"---- {n} overlap issues")


_save = Figure.savefig
_done = set()


def savefig(self, *a, **k):
    if id(self) not in _done:
        _done.add(id(self)); check(self)
    return _save(self, *a, **k)


Figure.savefig = savefig
if __name__ == "__main__":
    script = sys.argv[1]
    sys.argv = sys.argv[1:]
    runpy.run_path(script, run_name="__main__")
