


import xarray
# from xarray.core import indexing
from xarray.backends.common import BACKEND_ENTRYPOINTS, BackendEntrypoint #, _normalize_path, AbstractDataStore #, BackendArray
# from xarray.backends import CachingFileManager
from xarray.backends.locks import SerializableLock #, get_write_lock, ensure_lock
# from xarray.core.utils import FrozenDict
import h5py

#from . import psc

class RunInfo:
    def __init__(self, filename):
        self._file = h5py.File(filename)
        self.read_mfields()
        self.read_mcrds()
        #print(f'fields {self._fields}')
        #print(f'crds {self._crds}')
    
    def close(self):
        self._file.close()
        
    def fields(self):
        return self._fields
    
    def crds(self):
        return self._crds
    
    def find_pfx(self, pfx):
        for k in self._file:
            if k.startswith(pfx):
                return k

    def read_mfields(self):
        self._fields = {}
        for k in self._file:
            if k.startswith('jeh-') or k.startswith('all_1st-') or k.startswith('fields_vpic-') or k.startswith('hydro_vpic-'):
                mrc_fld = self._file[k]
                for kk in mrc_fld:
                    group = mrc_fld[kk]
                    dsets = [group[p]['3d'] for p in group]
                    self._fields[kk] = dsets
            if k.startswith('mrc_fld_'):
                mrc_fld = self._file[k]
                for kk in mrc_fld:
                    group = mrc_fld[kk]
                    dsets = [group[p]['3d'] for p in group]
                    self._fields[kk] = dsets
                    
    def read_mcrds(self, bnd=0):
        mcrds = {}
        for d, crd in enumerate(['x', 'y', 'z']):
            mcrds_path = self.find_pfx('crd[{}]-'.format(d))
            if mcrds_path is None:
                return self.read_mcrds_2(bnd)
            group = self._file[mcrds_path + '/crd[{}]'.format(d)]
            paths = [p + "/1d" for p in group]
            dsets = [group[path] for path in paths]
            if (bnd):
                dsets = [dset[bnd:-bnd] for dset in dsets]
            mcrds[crd] = dsets
        
        self._crds = mcrds                    
        
    def read_mcrds_2(self, bnd):
        bnd = 2 # FIXME
        mrc_crds = self._file[self.find_pfx('mrc_crds')]

        mcrds = {}
        for d, crd in enumerate(['x', 'y', 'z']):
            mcrds_path = 'mrc_m1_{}/crd{}'.format(mrc_crds.attrs['mcrd{}'.format(d)][0], d)
            group = self._file[mcrds_path]
            paths = [p + "/1d" for p in group]
            dsets = [group[path] for path in paths]
            if (bnd):
                dsets = [dset[bnd:-bnd] for dset in dsets]
            mcrds[crd] = dsets

        self._crds = mcrds
    
def psc_open_dataset(filename_or_obj, drop_variables=None):
    run = RunInfo(filename_or_obj)
    fields = run.fields()
    crds = run.crds()
    ds = xarray.Dataset({fieldname: xarray.DataArray(fields[fieldname][0][:].T,
                                                     coords=[crds[d][0][:] for d in 'xyz'],
                                                     dims=['x', 'y', 'z']) for fieldname in fields})
    ds.set_close(run.close)
    return ds

class PscHdf5BackendEntrypoint(BackendEntrypoint):
    available=True
    def open_dataset(
        self,
        filename_or_obj,
        *,
        drop_variables=None,
    ):
        return psc_open_dataset(filename_or_obj, drop_variables=drop_variables)

    open_dataset_parameters = ["filename_or_obj", "drop_variables"]

    def guess_can_open(self, filename_or_obj):
        try:
            _, ext = os.path.splitext(filename_or_obj)
        except TypeError:
            return False
        return ext in {".h5"}


if xarray.__version__ == "2023.4.1":
    # FIXME determine exactly when the API changed
    BACKEND_ENTRYPOINTS["pschdf5"] = ("psc", PscHdf5BackendEntrypoint)
else:
    # API of version 0.19.0
    BACKEND_ENTRYPOINTS["pschdf5"] = PscHdf5BackendEntrypoint
