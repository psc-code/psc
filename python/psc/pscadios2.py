
from .psc import Psc, FieldToComponent
from . import adios2py

import xarray as xr
from xarray.backends.common import BACKEND_ENTRYPOINTS
from xarray.backends import BackendEntrypoint, BackendArray
from xarray.core import indexing

from collections import namedtuple

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

class File:   
    def __init__(self, filename, length=None):
        # FIXME, opens and closes the file, ie, slows things unnecessarily
        self.psc = Psc(filename, length=length)
        
        self._fields_to_index = FieldToComponent(['he_e', 'e', 'i'])

        self._file = adios2py.file(filename)
        
        self._fields = {}
        for var in self._file.vars:
            for field, idx in self._fields_to_index[var].items():
                self._fields[field] = _FieldInfo(varname=var, component=idx)

    def __del__(self):
        pass
        #print("DBG: File closing")
        #self._file.close()

    @property
    def fields(self):
        return self._fields.keys()


def psc_open_dataset(filename, length=None, drop_variables=None):
    file = File(filename, length)
    fields = file.fields
    vars = {}
    for f in fields:
        field = file._fields[f]
        var = file._file[field.varname]
        coords = { "x": file.psc.x,
                   "y": file.psc.y,
                   "z": file.psc.z }
        arr = PscAdios2Array(var, field.component)
        vars[f] = xr.DataArray(arr, dims=['x', 'y', 'z'], coords=coords)

    return xr.Dataset(vars)

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

