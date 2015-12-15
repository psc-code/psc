#! /usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import re

def _find_path(file, obj_name):
    for path in file:
        expr = "^%s-uid-0x([0-9a-f]*)$" % (obj_name)
        match = re.match(expr, path)
        #print "path", path, "--", "%s-uid" % (obj_name), "expr", expr, match
        if match:
            return path
    raise Exception("object '%s' not found!" % obj_name)

class GgcmFields:
    def __init__(self, path, step, pfx="p"):
        filename = "%s/%s.3d.%06d_p%06d.h5" % (path, pfx, step, 0)
        print "Opening '%s'" % (filename)
        self._h5file = h5py.File(filename, 'r')

        self.crd = [self._read_crd(d) for d in xrange(3)]

    def _read_crd(self, dim):
        #path = "crd/%s-0" % ("xyz"[dim])
        #dset = self._h5file[path]
        path = _find_path(self._h5file, "crd\[%d\]" % dim)
        dset = self._h5file["%s/crd[%d]/p0/1d" % (path, dim)]
        crd = .5*(dset[1:-2] + dset[2:-1])

        return crd

    def _read_f3(self, field, comp):
        field = _find_path(self._h5file, field)
        #dset = self._h5file["%s/%s/3d" % (field, comp)]
        dset = self._h5file["%s/%s/p0/3d" % (field, comp)]
        fld = dset[:]

        return fld

    def __getitem__(self, what):
        if what in ["rr", "pp", "vx", "vy", "vz", "bx", "by", "bz"]:
            return self._read_f3("mrc_fld_" + what, what)
        elif what in ["rr1", "rv1x", "rv1y", "rv1z", "uu1", "b1x", "b1y", "b1z"]:
            return self._read_f3("mrc_fld_" + what, what)
        elif what in ["jx", "jy", "jz"]:
            return self._read_f3("mrc_fld_" + what, what)
        elif what in ["bdipx", "bdipy", "bdipz", "bdipccx", "bdipccy", "bdipccz"]:
            return self._read_f3("mrc_fld_" + what, what)
        elif what in ["ymask", "zmask", "rmask"]:
            return self._read_f3("mrc_fld_" + what, what)
        elif what in ["divB"]:
            return self._read_f3("mrc_fld", what)
        else:
            func = "_get_" + what
            return getattr(self, func)()

        raise KeyError

    def _get_bx(self):
        assert False
        b1x = self["b1x"]
        bx = np.empty_like(b1x)
        bx[:-1,:,:] = .5 * (b1x[:-1,:,:] + b1x[1:,:,:])
        return bx

    def _get_epar(self):
        hx, hy, hz = self["hx"], self["hy"], self["hz"]
        ex, ey, ez = self["ex"], self["ey"], self["ez"]
        return (ex * hx + ey * hy + ez * hz) * (hx**2 + hy**2 + hz**2)**-.5

    def _get_psi(flds):
        hz, hy = flds["hz"], flds["hy"]
        dz = flds.crd[2][1] - flds.crd[2][0]
        dy = flds.crd[1][1] - flds.crd[1][0]
        nz, ny, _ = hy.shape

        A = np.empty_like(hy).reshape(nz, ny)
        hz = hz.reshape(nz, ny)
        hy = hy.reshape(nz, ny)
        A[0,0] = 0.
        for i in range(1,nz):
            A[i,0] = A[i-1,0] + dz * ( hy[i,0] + hy[i-1,0] )/2.

        for j in range(1,ny):
            A[:,j] = A[:,j-1] - dy * ( hz[:,j-1] + hz[:,j] )/2.

        return A.reshape(nz, ny, 1)

def plot_slice_xy(flds, fld, xlim=None, ylim=None):
    crd = flds.crd
    if not xlim:
        xlim = (crd[0][0], crd[0][-1])
    if not ylim:
        ylim = (crd[1][0], crd[1][-1])
    mx, my, mz = fld.shape
    plt.pcolormesh(crd[0], crd[1], fld[mx/2-1,:,:])#, cmap="gray")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.gca().set_aspect("equal")
    plt.colorbar()

def plot_slice_xz(flds, fld, xlim=None, zlim=None):
    crd = flds.crd
    if not xlim:
        xlim = (crd[0][0], crd[0][-1])
    if not zlim:
        zlim = (crd[2][0], crd[2][-1])
    mx, my, mz = fld.shape
    plt.pcolormesh(crd[0], crd[2], fld[:,my/2-1,:])
    plt.xlim(xlim)
    plt.ylim(zlim)
    plt.gca().set_aspect("equal")
    plt.colorbar()

xlim = None
ylim = None
zlim = None

pfx = "run"
fldname = "rr"
step = 10
flds = GgcmFields(".", step, pfx)
fld = flds[fldname]# - GgcmFields(".", 0, pfx)[fldname]

plt.figure()
plot_slice_xy(flds, fld, xlim, zlim)
# circle = plt.Circle((0,0), 1, ec='k', fc='none')
# plt.gca().add_artist(circle)
#plt.clim(139.9, 140.1)
#plt.clim(-.2, .2)
plt.savefig("xy.png", dpi=150)
