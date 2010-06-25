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

def plot_z(ax, file, x, fld):
    ax.plot(x[:,0,0])

def plot_xz(ax, file, x, fld, jy=0):
    pcolormesh(x[:,jy,:])
    colorbar()
    xlim(0,x.shape[2])
    ylim(0,x.shape[0])
    
def main():
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] filename")
    parser.add_option("-f", "--fields", default="ne",
                      help="specify which fields to plot")
    parser.add_option("-d", "--diff",
                      help="diff to specified file")
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("specify exactly one filename")
    filename = args[0]
    file = h5py.File(filename, "r")
    if options.diff:
        file2 = h5py.File(options.diff, "r")

    flds = options.fields.split(",")
    n = len(flds)
    
    for i, fld in enumerate(flds):
        ax = subplot(n, 1, i+1)
        x = load_array(file['/psc/fields/%s' % fld])
        if options.diff:
            x2 = load_array(file2['/psc/fields/%s' % fld])
            x -= x2
            
        dims = x.shape[2], x.shape[1], x.shape[0]
        if dims[0] == 1 and dims[1] == 1:
            plot_z(ax, file, x, fld)
        elif dims[1] == 1:
            plot_xz(ax, file, x, fld)
        else:
            plot_xz(ax, file, x, fld, jy=x.shape[1]/2)
            #raise Exception("dims %s not supported" % (dims,))

    show()

if __name__ == "__main__":
    main()
    
