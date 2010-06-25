#! /usr/bin/env python

import os
import sys
import numpy
import numpy.linalg as linalg
import h5py

def load_array(ds):
    a = numpy.empty(shape=ds.shape, dtype=ds.dtype)
    a[:] = ds[:]
    return a

def main():
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] file1 file2")
    parser.add_option("-f", "--fields", default="ne",
                      help="specify which fields to plot")
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("specify exactly one filename")
    file1 = h5py.File(args[0], "r")
    file2 = h5py.File(args[1], "r")

    flds = options.fields.split(",")
    n = len(flds)
    
    for i, fld in enumerate(flds):
        x1 = load_array(file1['/psc/fields/%s' % fld])
        x2 = load_array(file2['/psc/fields/%s' % fld])
        diff = x1 - x2
        norm_inf = linalg.norm(diff.flatten(), numpy.inf)
        print "%s diff = %g" % (fld, norm_inf)

if __name__ == "__main__":
    main()
    
