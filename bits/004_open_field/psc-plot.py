#! /usr/bin/env python

import psc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

path = '.'

psc.downscale = 0

def plot_fld(fldname, step, cra=None, xlim=None, ylim=None):
    flds = psc.PscFields(path, step, "p")
    fld = flds[fldname]
    crd_nc = flds.crd_nc

    plt.figure()
    kwargs = {}

    plt.pcolormesh(crd_nc[1], crd_nc[2], fld[:,:,0], **kwargs)
    if cra:
        plt.clim(*cra)
    plt.colorbar()
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
    plt.savefig("%s-%06d.png" % (fldname, step), dpi=100)
    plt.close()


kwargs = {}
steps = xrange(0, 100, 1)

for step in steps:
    plot_fld("ex", step, cra=(-1, 1), **kwargs)
    plot_fld("ey", step, cra=(-1, 1), **kwargs)
    plot_fld("ex", step, cra=(-1, 1), **kwargs)
    plot_fld("hx", step, cra=(-1, 1), **kwargs)
    plot_fld("hy", step, cra=(-1, 1), **kwargs)
    plot_fld("hz", step, cra=(-1, 1), **kwargs)
