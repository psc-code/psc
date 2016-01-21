#! /usr/bin/env python

import psc
import numpy as np
import matplotlib.pyplot as plt

path = '.'

psc.downscale = 0
psc.d_i = 5.

def plot_fld(fldname, step, cra=None, plot_B=False, plot_log=False, xlim=None, ylim=None):
    flds = psc.PscFields(path, step, "p")
    fld = flds[fldname]
    crd = flds.crd

    if plot_log:
        fld = np.log10(fld)

    plt.figure()
    plt.pcolormesh(crd[2], crd[1], fld[:,:,0].transpose())
    if cra:
        plt.clim(*cra)
    plt.colorbar()
    if plot_B:
        psi = flds["psi"]
        plt.contour(crd[2], crd[1], psi[:,:,0].transpose(), 20, linestyles = 'solid', colors ='k')

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
steps = xrange(0, 100000, 1000)
#kwargs = { "xlim": (5, 10) }

for step in steps:
    plot_fld("jx", step, cra=[0,.15], plot_B=True, **kwargs)
