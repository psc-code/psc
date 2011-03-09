#! /usr/bin/env python

import h5py
import numpy
from scipy.io import savemat

def _make_filename(basename, time, proc):
    return "%s.%06d_p%06d.h5" % (basename, time, proc)

def _read_patches(obj, fld_name, domain, fld):
    obj_fld = obj[fld_name]
    nr_patches = obj.attrs['nr_patches']
    for p in xrange(nr_patches):
        obj_patch = obj_fld['p%d' % p]
        dset = numpy.transpose(obj_patch['3d'])
        shape = dset.shape
        gp = obj_patch.attrs['global_patch']
        off = domain['p%d' % gp].attrs['off']
        fld[off[0]:off[0]+shape[0],
            off[1]:off[1]+shape[1],
            off[2]:off[2]+shape[2]] = dset

def hdf5_to_mat(basename, times):
    print 'basename', basename
    print 'times', times

    flds = {}

    filename = _make_filename(basename, times[0], 0)
    h5 = h5py.File(filename, 'r')
    for fld_name in h5:
        obj = h5[fld_name]
        if obj.attrs['mrc_obj_class'] != 'mrc_m3':
            continue

        assert obj.attrs['nr_comps'] == 1
        domain_path = obj.attrs['domain']
        domain = h5[domain_path]
        gdims = domain.attrs['m']
        mpi_size = domain.attrs['mpi_size']
        if fld_name not in flds:
            flds[fld_name] = \
                numpy.zeros((len(times), gdims[0], gdims[1], gdims[2]))

    for time_idx, time in enumerate(times):
        for proc in xrange(0, mpi_size):
            filename = _make_filename(basename, time, proc)
            print 'filename', filename
            h5 = h5py.File(filename, 'r')
            for fld_name in flds.keys():
                obj = h5[fld_name]
                fld = flds[fld_name]
                # read each patch and put into global field at the right position
                _read_patches(obj, fld_name, domain, fld[time_idx])
            h5.close()
                
    for key, fld in flds.iteritems():
        print 'key', key, ':', fld.shape

    return flds

    # import pylab
    # pylab.pcolormesh(flds['jz'][0,:,:,0])
    # pylab.show()

if __name__ == "__main__":
    import sys
    flds = hdf5_to_mat(sys.argv[1], map(int, sys.argv[2:]))
    savemat("%s.mat" % sys.argv[1], flds)
