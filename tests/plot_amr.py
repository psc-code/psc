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
buf = 1
#patches = [0,1,2,4,7,11,12] #,13]
#patches = [5,6,9,10]
#patches = [13,14,16,17]
#patches = [5,6,9,10,13,14,16,17]
patches = [0,1,2,4,7,11,12,13, 5,6,9,10]
patches = xrange(22)
#patches = [0]
times = xrange(1,100,5)

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

def plot_component(basename, fldname, compname, time, symbol, **kwargs):
    for p in patches:
        fld, crdcc, crdnc = read_patch(basename, fldname, compname, time, p)
        if compname == "EX":
            slx = slicex(sw-buf, -(sw-buf))
            sly = slicex(sw-buf, -(sw-buf)+1)
            X, Y = np.meshgrid(crdcc[0][slx], crdnc[1][sly])
        elif compname == "EY":
            slx = slicex(sw-buf, -(sw-buf)+1)
            sly = slicex(sw-buf, -(sw-buf))
            X, Y = np.meshgrid(crdnc[0][slx], crdcc[1][sly])
        elif compname == "EZ":
            slx = slicex(sw-buf, -(sw-buf)+1)
            sly = slicex(sw-buf, -(sw-buf)+1)
            X, Y = np.meshgrid(crdnc[0][slx], crdnc[1][sly])
        elif compname == "HY":
            slx = slicex(sw-buf, -(sw-buf))
            sly = slicex(sw-buf, -(sw-buf)+1)
            X, Y = np.meshgrid(crdcc[0][slx], crdnc[1][sly])
        elif compname == "HZ":
            slx = slicex(sw-buf, -(sw-buf))
            sly = slicex(sw-buf, -(sw-buf))
            X, Y = np.meshgrid(crdcc[0][slx], crdcc[1][sly])

        fld = fld[sw,sly,slx]

        ax = plt.gca(projection='3d')
        if '-' in symbol:
            ax.plot_wireframe(X, Y, fld, **kwargs)
        if '.' in symbol:
            ax.plot(X.flatten(), Y.flatten(), fld.flatten(), '.', **kwargs)
        #ax.plot(X[2,:], fld[2,:,m], 'o-', **kwargs) # 2d
        # if buf > 0:
        #     ax.plot_wireframe(X[buf:-buf,buf:-buf], Y[buf:-buf,buf:-buf], fld[buf:-buf,buf:-buf,m], color='g')
        #ax.plot_wireframe(X, Y, np.sin(.5+2*np.pi*X)*np.cos(.5+2*np.pi*Y), color='g')

def movie():
    global buf
    buf = 0
    plt.ion()
    plt.figure(figsize=(16,10))
    for time in times:
        plt.clf()
        plot_component(basename, "fld", "EX", time, '-', color='r')
        plt.draw()
        #plt.show()
        
    plt.show()

def boundary():
    #plot_component(basename, "fld", "EZ", 0, '.', color='r')
    plot_component(basename, "fld", "EX", 1, '-', color='r')
    plt.show()

boundary()
#movie()
