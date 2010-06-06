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

def plot_z(ax, filename, fld):
    f = h5py.File(filename, "r")
    x = load_array(f['/psc/fields/%s' % fld])
    # X = -load_array(f['mhd/crd/x'])
    # Y = load_array(f['mhd/crd/y'])
    # Z = load_array(f['mhd/crd/z'])

    print x.shape
    plot(x[:,0,-0])

def main():
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] filename")
    parser.add_option("-f", "--fields", default="ne",
                      help="specify which fields to plot")
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("specify exactly one filename")
    filename = args[0]

    flds = options.fields.split(",")
    n = len(flds)
    for i, fld in enumerate(flds):
        ax = subplot(n, 1, i+1)
        plot_z(ax, filename, fld)

    show()

if __name__ == "__main__":
    main()
    
