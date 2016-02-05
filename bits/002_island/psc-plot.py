#! /usr/bin/env python

import psc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

path = '.'

psc.downscale = 0
psc.d_i = 5.

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def plot_fld(fldname, step, cra=None, plot_B=False, plot_log=False, xlim=None, ylim=None, midpoint=None, vmin=None, vmax=None):
    flds = psc.PscFields(path, step, "p")
    fld = flds[fldname]
    crd = flds.crd

    if plot_log:
        fld = np.log10(fld)

    plt.figure()
    kwargs = {}
    if midpoint is not None:
        kwargs = { "norm": MidpointNormalize(vmin=vmin, vmax=vmax, midpoint=midpoint),"cmap": plt.cm.seismic }

    plt.pcolormesh(crd[2], crd[1], fld[:,:,0].transpose(), **kwargs)
    if cra:
        plt.clim(*cra)
    plt.colorbar()
    if plot_B:
        psi = flds["psi"]
        plt.contour(crd[2], crd[1], psi[:,:,0].transpose(), 10, linestyles = 'dashed', colors ='k')

    if xlim:
        plt.xlim(*xlim)
    else:
        plt.xlim(crd[2][0], crd[2][-1])

    if ylim:
        plt.ylim(*ylim)
    else:
        plt.ylim(crd[1][0], crd[1][-1])

    plt.axes().set_aspect('equal')
    plt.savefig("%s-%06d.png" % (fldname, step), dpi=300)
    plt.close()


kwargs = {}
steps = xrange(00000, 100000, 1000)
#kwargs = { "xlim": (5, 10) }

for step in steps:
    plot_fld("jx", step, plot_B=True, midpoint=0, vmin=-0.3, vmax = .1, **kwargs)
    # plot_fld("hz", step, **kwargs)
    # plot_fld("vx_e", step, **kwargs)
    # plot_fld("vx_i", step, **kwargs)
    #    plot_fld("divb", step, **kwargs)
