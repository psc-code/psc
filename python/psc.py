#! /usr/bin/env python

import h5py
import numpy
from pylab import *

downscale = 0
d_i = 1 # used to rescale from d_e (natural normalization) to d_i

def _find_path(file, obj_name):
    for path in file:
        if path.startswith("%s-uid" % (obj_name)):
            return path
    raise Exception("object '%s' not found!" % obj_name)

class PscFields:
    def __init__(self, path, step, pfx="p"):
        filename = "%s/%sfd.%06d_p%06d.h5" % (path, pfx, step, 0)
        print "Opening '%s'" % (filename)
        self._h5file = h5py.File(filename, 'r')

        self.crd = [self._read_crd(d) / d_i for d in xrange(3)]

        path = _find_path(self._h5file, "psc")
        self.cc = self._h5file[path].attrs["cc"][0]
        self.time = self._h5file[path].attrs["time"][0] / self.cc
        self.timestep = self._h5file[path].attrs["timestep"][0]

    def _read_crd(self, dim):
        path = _find_path(self._h5file, "crd[%d]" % (dim))
        dset = self._h5file["%s/crd[%d]/p0/1d" % (path, dim)]
        sw = self._h5file[path].attrs["sw"][0]
        crd = dset[:]

        for i in xrange(downscale):
            crd = .5*(crd[::2] + crd[1::2])

        return crd

    def _read_f3(self, field, comp):
        field = _find_path(self._h5file, field)
        dset = self._h5file["%s/%s/p0/3d" % (field, comp)]
        fld = dset[:]

        for i in xrange(downscale):
            fld = .25*(fld[::2,::2] + fld[1::2,::2] + fld[::2,1::2] + fld[1::2,1::2])

        return fld

    def __getitem__(self, what):
        if what.startswith("n_"):
            try:
                return self._read_f3("n_1st_single", what)
            except:
                return self._read_f3("n_1st_c", what)
        if what.startswith("vx_") or what.startswith("vy_") or what.startswith("vz_"):
            try:
                return self._read_f3("v_1st_single", what)
            except:
                return self._read_f3("v_1st_c", what)
        elif what in ["ex", "ey", "ez"]:
            return self._read_f3("e", what)
        elif what in ["hx", "hy", "hz"]:
            return self._read_f3("h", what)
        elif what in ["jx", "jy", "jz"]:
            return self._read_f3("j", what)
        elif what in ["ne", "ni", "nn"]:
            return self._read_f3("n", what)
        elif what in ["dive", "divj", "divb"]:
            return self._read_f3(what, what)
        else:
            func = "_get_" + what
            return getattr(self, func)()

        raise KeyError

    def _get_epar(self):
        hx, hy, hz = self["hx"], self["hy"], self["hz"]
        ex, ey, ez = self["ex"], self["ey"], self["ez"]
        return (ex * hx + ey * hy + ez * hz) * (hx**2 + hy**2 + hz**2)**-.5

    def _get_h_tot(self):
        hx, hy, hz = self["hx"], self["hy"], self["hz"]
        return (hx**2 + hy**2 + hz**2) ** .5
        
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

        
class PscParticles:
    def __init__(self, path, step):
        filename = "%s/prt.%06d_p%06d.h5" % (path, step, 0)
        print "Opening '%s'" % (filename)
        self._h5file = h5py.File(filename, 'r')

        # path = _find_path(self._h5file, "psc")
        # self.time = self._h5file[path].attrs["time"][0]
        # self.timestep = self._h5file[path].attrs["timestep"][0]

        self.data = self._h5file["particles/p0/1d"]

