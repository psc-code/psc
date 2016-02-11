#! /usr/bin/env python

import psc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

path = '.'

psc.downscale = 0

def plot_fld(fldname, step, cra=None, plot_B=False, xlim=None, ylim=None):
    flds = psc.PscFields(path, step, "p")
    fld = flds[fldname]
    crd_nc = flds.crd_nc
    crd_cc = flds.crd

    plt.figure()
    kwargs = {}

    plt.pcolormesh(crd_nc[1], crd_nc[2], fld[:,:,0], **kwargs)
    if cra:
        plt.clim(*cra)
    plt.colorbar()

    if plot_B:
        psi = flds["psi"]
        plt.contour(crd_cc[1], crd_cc[2], psi[:,:,0], 20, linestyles = 'solid', colors ='k')

    plt.title(fldname)
    plt.xlabel("y")
    plt.ylabel("z")

    if xlim:
        plt.xlim(*xlim)
    else:
        plt.xlim(crd_nc[1][0], crd_nc[1][-1])

    if ylim:
        plt.ylim(*ylim)
    else:
        plt.ylim(crd_nc[2][0], crd_nc[2][-1])

    plt.axes().set_aspect('equal')
    plt.savefig("%s-%06d.png" % (fldname, step), dpi=300)
    plt.close()


kwargs = {}
steps = xrange(000, 1000, 10)

for step in steps:
    plot_fld("n_e", step, cra=[.7, 1.3], plot_B=False, **kwargs)
    plot_fld("n_i", step, cra=[.7, 1.3], plot_B=False, **kwargs)
    # plot_fld("hx", step, **kwargs)
    # plot_fld("hy", step, **kwargs)
    # plot_fld("hz", step, **kwargs)
    # plot_fld("jx", step, **kwargs)
    # plot_fld("jy", step, **kwargs)
    # plot_fld("jz", step, **kwargs)
