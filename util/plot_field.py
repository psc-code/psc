#! /usr/bin/env python

import os
import sys
import numpy
import h5py
from pylab import *

def load_array(ds):
    a = numpy.empty(shape=ds.shape, dtype=ds.dtype)
    a[:] = ds[:]
    return a

def plot_z(ax, file, fld):
    x = load_array(file['/psc/fields/%s' % fld])

    ax.plot(x[:,0,0])

def plot_xz(ax, file, fld):
    x = load_array(file['/psc/fields/%s' % fld])

    ax.pcolormesh(x[:,0,:])

def main():
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] filename")
    parser.add_option("-f", "--fields", default="ne",
                      help="specify which fields to plot")
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("specify exactly one filename")
    filename = args[0]
    file = h5py.File(filename, "r")
    lo = file["psc"].attrs["lo"]
    hi = file["psc"].attrs["hi"]
    dims = hi - lo

    flds = options.fields.split(",")
    n = len(flds)
    
    for i, fld in enumerate(flds):
        ax = subplot(n, 1, i+1)
        if dims[0] == 1 and dims[1] == 1:
            plot_z(ax, file, fld)
        elif dims[1] == 1:
            plot_xz(ax, file, fld)
        else:
            raise Exception("dims %s not supported" % dims)

    show()

if __name__ == "__main__":
    main()
    
