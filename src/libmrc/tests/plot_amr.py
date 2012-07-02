#! /usr/bin/env python

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

mx = 8
my = 8
nr_patches = 10

def read_patch(basename, fldname, time, p):
    filename = "%s.%06d_p%06d_%s.asc" % (basename, time, p, fldname)
    fld = np.loadtxt(filename)
    assert fld.shape[0] == mx * my
    nr_comps = fld.shape[1] - 3
    crdcc = [fld[:mx,0], fld[::mx,1]]
    crdnc = []
    for crd in crdcc:
        tmpcc = np.empty(crd.size + 2)
        tmpcc[1:-1] = crd
        tmpcc[0] = 2 * tmpcc[1] - tmpcc[2]
        tmpcc[-1] = 2 * tmpcc[-2] - tmpcc[-3]
        crdnc.append(.5 * (tmpcc[:-2] + tmpcc[1:-1]))
    fld = fld[:,3:].reshape(mx, my, nr_comps)
    return fld, crdcc, crdnc

def plot_component(basename, fldname, time, m, **kwargs):
    for p in xrange(nr_patches):
        fld, crdcc, crdnc = read_patch(basename, fldname, time, p)
        if m == 1: # EY
            X, Y = np.meshgrid(crdnc[0], crdcc[1])
        elif m == 5: # BZ
            X, Y = np.meshgrid(crdcc[0], crdcc[1])

        ax = plt.gca(projection='3d')
        ax.plot_wireframe(X, Y, fld[:,:,m], **kwargs)

for time in xrange(2):
    plt.figure()
    plot_component("run", "fld", time, 1, color='r')
    plot_component("run", "fld", time, 5, color='b')

plt.show()


