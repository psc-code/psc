#! /usr/bin/env python

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import h5py

basename = "run"
mx = 8
my = 8
mz = 1
sw = 3
buf = 0
#patches = [0,1,2,4,7,11,12,13]
#patches = [5,6,9,10]
#patches = [0,1,2,4,7,11,12,13, 5,6,9,10]
patches = xrange(10)
times = xrange(0,32,1)

EX = 0
EY = 1
EZ = 2
HX = 3
HY = 4
HZ = 5

def read_patch(basename, fldname, compname, time, p):
    rank = 0
    filename = "%s.%06d_p%06d.h5" % (basename, time, rank)
    f = h5py.File(filename, 'r')
    dset = f['%s/%s/p%d/3d' % (fldname, compname, p)]
    fld = dset[:,:,:]
    assert fld.shape[2] == mx + 2 * sw
    assert fld.shape[1] == my + 2 * sw
    assert fld.shape[0] == mz + 2 * sw
    crdnc = [f['crd0/p%d/1d' % p][:], f['crd1/p%d/1d' % p][:]]
    crdcc = [.5*(c[1:] + c[:-1]) for c in crdnc]
    crdnc = [c[:-1] for c in crdnc]
    return fld, crdcc, crdnc

def slicex(l, h):
    if h == 0:
        return slice(l, None)
    else:
        return slice(l, h)

def plot_component(basename, fldname, compname, time, **kwargs):
    for p in patches:
        fld, crdcc, crdnc = read_patch(basename, fldname, compname, time, p)
        if compname == "EY":
            slx = slicex(sw-buf, -(sw-buf)+1)
            sly = slicex(sw-buf, -(sw-buf))
            X, Y = np.meshgrid(crdnc[0][slx], crdcc[1][sly])
        elif compname == "EZ":
            slx = slicex(sw-buf, -(sw-buf)+1)
            sly = slicex(sw-buf, -(sw-buf)+1)
            X, Y = np.meshgrid(crdnc[0][slx], crdnc[1][sly])
        elif compname == "HZ":
            slx = slicex(sw-buf, -(sw-buf))
            sly = slicex(sw-buf, -(sw-buf))
            X, Y = np.meshgrid(crdcc[0][slx], crdcc[1][slx])

        fld = fld[sw,sly,slx]

        ax = plt.gca(projection='3d')
        ax.plot_wireframe(X, Y, fld, **kwargs)
        #ax.plot(X[2,:], fld[2,:,m], 'o-', **kwargs)
        # if buf > 0:
        #     ax.plot_wireframe(X[buf:-buf,buf:-buf], Y[buf:-buf,buf:-buf], fld[buf:-buf,buf:-buf,m], color='g')
        #ax.plot_wireframe(X, Y, np.sin(.5+2*np.pi*X)*np.cos(.5+2*np.pi*Y), color='g')

plt.ion()
for time in times:
    #plt.figure()
    plt.clf()
    plot_component(basename, "fld", "EY", time, color='r')
    #plot_component("run", "fld", "EZ", time, color='r')
    #plot_component(basename, "fld", "HZ", time, color='b')
    plt.draw()
    #plt.show()

plt.show()


