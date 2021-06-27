
from .psc import RunInfo, FieldToComponent
from . import adios2py

import xarray
from xarray.core import indexing
from xarray.backends.common import BACKEND_ENTRYPOINTS, BackendEntrypoint, BackendArray, _normalize_path, AbstractDataStore
from xarray.backends import CachingFileManager

from collections import namedtuple
import logging

class PscAdios2Array(BackendArray):
    """Lazy evaluation of a variable stored in PSC's adios2 field output.
    
    This takes care of slicing out the specific component of the data stored as 4-d array.
    """
    def __init__(self, var, m):
        self._var = var
        self._m = m
        self.shape = var.shape[:-1]
        self.dtype = var.dtype
        
    def __getitem__(self, key):
        return indexing.explicit_indexing_adapter(
            key, self.shape, indexing.IndexingSupport.BASIC, self._getitem)
    
    def _getitem(self, args):
        return self._var[(*args, self._m)]

_FieldInfo = namedtuple('FieldInfo', ['varname', 'component'])

class PscAdios2Store(AbstractDataStore):   
    def __init__(self, manager, length=None):
        self._manager = manager
        self.psc = RunInfo(self.ds, length=length)        
        
    @classmethod
    def open(cls, filename, length=None):
        manager = CachingFileManager(adios2py.file, filename)
        #manager = adios2py.file(filename)
        return cls(manager, length=length)
    
    def _acquire(self, needs_lock=True):
        with self._manager.acquire_context(needs_lock) as root:
            ds = root
        return ds

    @property
    def ds(self):
        return self._acquire()

    def get_variables(self):
        fields_to_index = FieldToComponent(['he_e', 'e', 'i'])

        vars = {}
        for var in self.ds.variables:
            for field, idx in fields_to_index[var].items():
                coords = { "x": self.psc.x, "y": self.psc.y, "z": self.psc.z }
                arr = PscAdios2Array(self.ds[var], idx)
                vars[field] = xarray.DataArray(arr, dims=['x', 'y', 'z'], coords=coords)

        return vars

    def get_attrs(self):
        return {}


def psc_open_dataset(filename_or_obj, length=None, drop_variables=None):
    filename_or_obj = _normalize_path(filename_or_obj)
    store = PscAdios2Store.open(filename_or_obj, length)
    
    vars, attrs = store.load()
    ds = xarray.Dataset(vars, attrs=attrs)
    ds.set_close(store.close)
    return ds

class PscAdios2BackendEntrypoint(BackendEntrypoint):
    def open_dataset(
        self,
        filename_or_obj,
        *,
        drop_variables=None,
        length=None
    ):
        return psc_open_dataset(filename_or_obj, drop_variables=drop_variables, length=length)

    open_dataset_parameters = ["filename_or_obj", "drop_variables"]

    def guess_can_open(self, filename_or_obj):
        try:
            _, ext = os.path.splitext(filename_or_obj)
        except TypeError:
            return False
        return ext in {".bp"}
    
BACKEND_ENTRYPOINTS["pscadios2"] = PscAdios2BackendEntrypoint

