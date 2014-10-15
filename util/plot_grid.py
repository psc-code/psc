#!/usr/bin/env python
"""Plot grid from grid2 files (RUN.grid2 and RUN.hgrid2)

Some --help is available for this script.

The grid can be plotted over data using the --sample -p and -o options.
In the plots, dashed lines represent [dxmin, 2.0 * dxmin, 10.0 * dxmin].
In addition, circles are plotted with radii [1.0, 3.0, 6.0, 10.0, 20.0].
"""

from __future__ import division, print_function
import sys
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# the matplotlib masked array stuff is too chatty
warnings.simplefilter("ignore")

# read grid
def read_grid(fname):
    crds = [None] * 3
    with open(fname, 'r') as f:
        for i in range(3):
            nx = int(f.readline())
            crdx = np.zeros(nx)
            for j in range(nx):
                crdx[j] = float(f.readline())
            crds[i] = crdx
    return crds

def plot_2d_gridcells(x, y, circles=None, colors='k', ax=None):
    if ax is None:
        ax = plt.gca()
    if circles is None:
        circles = []
    X, Y = np.meshgrid(x, y)
    blank = 0.0 * X + 0.0 * Y
    cmap = plt.cm.jet
    cmap.set_bad('w', 0.0)
    p = ax.pcolormesh(X, Y, np.ma.masked_where(blank == 0.0, blank),
                       edgecolors='k', linewidths=0.2,
                       antialiased=True, cmap=cmap)
    for i, radius in enumerate(circles):
        circle = plt.Circle((0, 0), radius, color=colors[i % len(colors)],
                            fill=False)
        ax.add_artist(circle)
    ax.axis('equal')
    return p

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s', "--slices", default="y=0.0,x=0.0",
                        help="comma separated list of slices, one figure per "
                             "slice; default is 'y=0.0,x=0.0'")
    parser.add_argument("--sample", default="",
                        help="optional datafile to plot sample data with grid "
                             "(requires Viscid to be in your PYTHONPATH)")
    parser.add_argument('-p', "--sample_var", default="rr",
                        help="optional sample quantity to plot")
    parser.add_argument('-o', "--plot_opts", default="",
                        help="options for sample plot")
    parser.add_argument("grid2", nargs=1, help="RUN.grid2 file")
    parser.add_argument("hgrid2", nargs=1, help="RUN.hgrid2 file")
    args = parser.parse_args()

    x, y, z = read_grid(args.grid2[0])
    dx, dy, dz = read_grid(args.hgrid2[0])

    dxmin = np.min(dx)
    dymin = np.min(dy)
    dzmin = np.min(dz)

    # setup fld if we're plotting a sample
    fld = None
    if args.sample:
        try:
            import viscid
            from viscid.plot import mpl
            viscid.readers.openggcm.GGCMGrid.mhd_to_gse_on_read = False
            fld = viscid.load_file(args.sample)[args.sample_var]
        except ImportError:
            print("Must have Viscid in your PYTHONPATH to plot a sample")
        except KeyError:
            print("Warning; unknown variable for sample")
        except TypeError:
            # given if the file doesnt exist, it will print its own warning
            pass
        except RuntimeError:
            # given if the calculator fails, it will print its own warning
            pass

    wrats = [1, 5]
    hrats = [1, 3]
    guides = [1.0, 2.0, 10.0]
    guide_colors = "rgy"
    circles = [1.0, 3.0, 6.0, 10.0, 20.0]
    circle_colors = 'krgby'

    for slc in args.slices.split(','):
        plane_dir, plane_loc = slc.lower().split('=')
        xstr, ystr = "xyz".replace(plane_dir, '')
        plane_str = "{0} Plane ({1} = {2})".format("-".join([xstr, ystr]).upper(),
                                                   plane_dir, plane_loc)

        _ = plt.figure()
        gspec = GridSpec(2, 2, width_ratios=wrats, height_ratios=hrats)
        ax_mesh = plt.subplot(gspec[3])
        ax_xcrds = plt.subplot(gspec[1], sharex=ax_mesh)
        ax_ycrds = plt.subplot(gspec[2], sharey=ax_mesh)

        # plot grid cells
        if fld is not None:
            mpl.plot(fld, slc, ax=ax_mesh, style="contourf", levels=100,
                     plot_opts=args.plot_opts,
                     colorbar=dict(ax=[ax_mesh, ax_xcrds], fraction=0.1))

        xarr = [x, y, z]["xyz".index(xstr)]
        yarr = [x, y, z]["xyz".index(ystr)]
        dxarr = [dx, dy, dz]["xyz".index(xstr)]
        dyarr = [dx, dy, dz]["xyz".index(ystr)]
        _dxmin = np.min(dxarr)
        _dymin = np.min(dyarr)
        _dxmax = np.max(dxarr)
        _dymax = np.max(dyarr)

        plot_2d_gridcells(xarr, yarr, circles=circles, colors=circle_colors,
                          ax=ax_mesh)
        ax_mesh.set_xlabel(xstr)

        # plot x vs dx for horizontal axis
        ax_xcrds.plot(xarr, dxarr)
        for i, guide in enumerate(guides):
            if guide * _dxmin > _dxmax:
                continue
            ax_xcrds.axhline(guide * _dxmin,
                             color=guide_colors[i % len(guide_colors)],
                             linestyle='--')
        ax_xcrds.set_ylabel('d' + xstr)

        # plot x vs dx for vertical axis
        ax_ycrds.plot(dyarr, yarr)
        for i, guide in enumerate(guides):
            if guide * _dymin > _dymax:
                continue
            ax_ycrds.axvline(guide * _dymin,
                             color=guide_colors[i % len(guide_colors)],
                             linestyle='--')
        ax_ycrds.set_ylabel(ystr)
        ax_ycrds.set_xlabel("d" + ystr)

        plt.suptitle(plane_str)
        info = """min dx: {0:.3g}
min dy: {1:.3g}
min dz: {2:.3g}""".format(dxmin, dymin, dzmin)
        ax = plt.subplot(gspec[0])
        ax.axis('off')
        ax.annotate(info, xy=(0, 0), xytext=(-0.3, 0.3),
                    textcoords='axes fraction')

    plt.show()

    return 0

if __name__ == "__main__":
    sys.exit(main())

##
## EOF
##
