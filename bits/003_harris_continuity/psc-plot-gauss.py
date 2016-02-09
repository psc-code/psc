#! /usr/bin/env python

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.size'] = 8
import psc
import numpy as np
import matplotlib.pyplot as plt

path = '.'

psc.downscale = 0
psc.d_i = 5.

def plot_one(ax, title, crd_nc, fld):
    ax.set_title(title)
    mappable = ax.pcolormesh(crd_nc[1], crd_nc[2], fld[:,:,0])
    ax.set_aspect('equal')
    ax.set_xlim(crd_nc[1][0], crd_nc[1][-1])
    ax.set_ylim(crd_nc[2][0], crd_nc[2][-1])
    cbar = plt.colorbar(mappable, ax=ax)

for step in xrange(0, 40, 10):
    flds = psc.PscFields(path, step, "continuity")
    crd_nc = flds.crd_nc
    div_j = flds["div_j"]
    d_rho = flds["d_rho"]

    f, (ax0, ax1, ax2) = plt.subplots(1, 3)

    plot_one(ax0, "d_rho", crd_nc, d_rho)
    plot_one(ax1, "div_j", crd_nc, div_j)
    plot_one(ax2, "sum"  , crd_nc, d_rho + div_j)

    plt.savefig("continuity-%06d.png" % (step), dpi=300)
    plt.close()

for step in xrange(0, 50, 10):
    flds = psc.PscFields(path, step, "gauss")
    crd_nc = flds.crd_nc
    rho = flds["rho"]
    dive = flds["div_E"]

    f, (ax0, ax1, ax2) = plt.subplots(1, 3)

    plot_one(ax0, "rho" , crd_nc, rho)
    plot_one(ax1, "dive", crd_nc, dive)
    diff = dive - rho
    # at conducting wall, Gauss's Law isn't expected to be satisfied,
    # so zero out the difference on the left (the right boundary isn't
    # included here, anyway)
    diff[:,0,:] = 0.
    print "gauss: max diff", np.max(np.abs(diff))
    plot_one(ax2, "diff", crd_nc, diff)

    plt.savefig("gauss-%06d.png" % (step), dpi=300)
    plt.close()
