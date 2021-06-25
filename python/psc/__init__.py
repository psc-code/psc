
import numpy as np
import psc.adios2py
import os
import xarray as xr
from collections import namedtuple

field_entry = namedtuple('field_entry', ['filename', 'varname', 'index'])

class Psc:
    def __init__(self, filename, length=None):
        file = adios2py.file(filename)
        assert len(file.vars) > 0
        var = next(iter(file.vars))
        self.gdims = np.asarray(file[var].shape)[0:3]
        file.close() # FIXME, should support with ... as

        if length is not None:
            self.length = np.asarray(length)
            # FIXME, corner should also be passed (or rather read)
            self.corner = -.5 * self.length
        else:
            self.length = self.gdims
            self.corner = np.array([0., 0., 0.])
            
        self.x = np.linspace(self.corner[0], self.corner[0] + self.length[0], self.gdims[0], endpoint=False)
        self.y = np.linspace(self.corner[1], self.corner[1] + self.length[1], self.gdims[1], endpoint=False)
        self.z = np.linspace(self.corner[2], self.corner[2] + self.length[2], self.gdims[2], endpoint=False)

    def __repr__(self):
        return f"Psc(gdims={self.gdims}, length={self.length}, corner={self.corner})"

class File:
    def __init__(self, filename, length=None):
        # FIXME, opens and closes the file, ie, slows things unnecessarily
        self.psc = Psc(filename, length=length)
        
        self._fields_to_index = FieldsToIndex(['he_e', 'e', 'i'])

        self._file = adios2py.file(filename)
        
        self._fields = {}
        for var in self._file.vars:
            for field, idx in self._fields_to_index[var].items():
                self._fields[field] = { 'var' : var, 'index' : idx }

    def __del__(self):
        print("DBG: File closing")
        self._file.close()

    @property
    def fields(self):
        return self._fields.keys()

    def read(self, fldname, start, count):
        field = self._fields[fldname]
        var = self._file[field['var']]
        arr = var[start[0]:start[0]+count[0],
                  start[1]:start[1]+count[1],
                  start[2]:start[2]+count[2],
                  field['index']]
        coords = { "x": self.psc.x[start[0]:start[0]+count[0]],
                   "y": self.psc.y[start[1]:start[1]+count[1]],
                   "z": self.psc.z[start[2]:start[2]+count[2]], }
        return xr.DataArray(arr, dims=['x', 'y', 'z'], coords=coords)


class FieldsToIndex:
    def __init__(self, species):
        self._map = {}
        self._map['jeh'] = { 'jx_ec': 0, 'jy_ec': 1, 'jz_ec' : 2,
                                         'ex_ec': 3, 'ey_ec': 4, 'ez_ec' : 5,
                                         'hx_fc': 6, 'hy_fc': 7, 'hz_fc' : 8 }
        self._map["all_1st"] = self._make_all_1st_to_index(species)
        
    def __getitem__(self, field):
        return self._map[field]

    def _make_all_1st_to_index(self, species):
        to_index = {}
        for i, s in enumerate(species):
            to_index[f'rho_{s}'] = 0 + 13 * i
            to_index[f'jx_{s}']  = 1 + 13 * i
            to_index[f'jy_{s}']  = 2 + 13 * i
            to_index[f'jz_{s}']  = 3 + 13 * i
            to_index[f'px_{s}']  = 4 + 13 * i
            to_index[f'py_{s}']  = 5 + 13 * i
            to_index[f'pz_{s}']  = 6 + 13 * i
            to_index[f'txx_{s}'] = 7 + 13 * i
            to_index[f'txx_{s}'] = 8 + 13 * i
            to_index[f'txx_{s}'] = 9 + 13 * i
            
        return to_index
    
        
class run:
    def __init__(self, path, L=None, pfx='pfd'):
        self._path = path
        self._pfx = pfx

        # open first var in first file to figure out global dims
        self._files = [f for f in os.listdir(path) if f.startswith(self._pfx) and f.endswith('.bp')]
        if not self._files:
            raise RuntimeError(f'No "{self._pfx}" data files found in "{path}"')
            
        self._files.sort()
        self._steps = {}
        for file in self._files:
            _, step, _ = file.split('.')
            step = int(step)
            if step in self._steps:
                self._steps[step].append(file)
            else:
                self._steps[step] = [file]
        # print("steps", self._steps.keys())

        self.psc = Psc(os.path.join(path, next(iter(self._files))), length=L)

        self._fields_to_index = FieldsToIndex(['he_e', 'e', 'i'])

        first_step = next(iter(self._steps))
        self.activate_time(first_step)
        
    def _make_all_1st_to_index(self, species):
        to_index = {}
        for i, s in enumerate(species):
            to_index[f'rho_{s}'] = 0 + 13 * i
            to_index[f'jx_{s}']  = 1 + 13 * i
            to_index[f'jy_{s}']  = 2 + 13 * i
            to_index[f'jz_{s}']  = 3 + 13 * i
            to_index[f'px_{s}']  = 4 + 13 * i
            to_index[f'py_{s}']  = 5 + 13 * i
            to_index[f'pz_{s}']  = 6 + 13 * i
            to_index[f'txx_{s}'] = 7 + 13 * i
            to_index[f'txx_{s}'] = 8 + 13 * i
            to_index[f'txx_{s}'] = 9 + 13 * i
            
        return to_index
    
    def activate_time(self, step):
        if not step in self._steps:
            raise ValueError(f'step {step} not found!')
        self._step = step
        filenames = self._steps[step]

        self.fields = {}
        for filename in filenames:
            file = adios2py.file(os.path.join(self._path, filename))
            for varname in file.vars:
                fields_to_index = self._fields_to_index[varname]
                # assert # of comps (and gdims?) match
                for f in fields_to_index:
                    self.fields[f] = field_entry(filename=filename, varname=varname, index=fields_to_index[f])
            file.close() # FIXME, should support with ... as 
        # print(f'fields {self.fields.keys()}')
        
    def read(self, fldname, start, count):
        field = self.fields[fldname]
        file = adios2py.file(os.path.join(self._path, field.filename))
        var = file[field.varname]
        m = field.index
        arr = var[start[0]:start[0]+count[0],
                  start[1]:start[1]+count[1],
                  start[2]:start[2]+count[2],
                  m]
        file.close()
        coords = { "x": self.psc.x[start[0]:start[0]+count[0]],
                   "y": self.psc.y[start[1]:start[1]+count[1]],
                   "z": self.psc.z[start[2]:start[2]+count[2]], }
        return xr.DataArray(arr, dims=['x', 'y', 'z'], coords=coords)

def read(run, time):
    py_e = run.read('py_he_e', [0, 0, 0], run.psc.gdims)
    n_e = -run.read('rho_he_e', [0, 0, 0], run.psc.gdims)
    return {"n_e" : n_e, "py_e" : py_e}