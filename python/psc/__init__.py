
import numpy as np
import psc.adios2py
import os
import xarray as xr
from collections import OrderedDict, namedtuple

field_entry = namedtuple('field_entry', ['filename', 'varname', 'index'])

class run:
    _jeh_to_index = { 'jx_ec': 0, 'jy_ec': 1, 'jz_ec' : 2,
                      'ex_ec': 3, 'ey_ec': 4, 'ez_ec' : 5,
                      'hx_fc': 6, 'hy_fc': 7, 'hz_fc' : 8 }

    _all_1st_to_index = { 'rho_he_e': 0, 'rho_e': 13, 'rho_i' : 26,
                          'px_he_e': 4, 'px_e': 17, 'px_i' : 30,
                          'py_he_e': 5, 'py_e': 18, 'py_i' : 31,
                          'txx_he_e': 7, 'txx_e': 20, 'txx_i' : 33,
                          'tyy_he_e': 8, 'tyy_e': 21, 'tyy_i' : 34,
                          'tzz_he_e': 9, 'tzz_e': 22, 'tzz_i' : 35,
                        }
    
    _fields_to_index = { "jeh": _jeh_to_index, "all_1st": _all_1st_to_index }
    
    def __init__(self, path, L=None, pfx='pfd'):
        self._path = path
        self._pfx = pfx

        # open first var in first file to figure out global dims
        self._files = [f for f in os.listdir(path) if f.startswith(self._pfx) and f.endswith('.bp')]
        if not self._files:
            raise RuntimeError(f'No "{self._pfx}" data files found in "{path}"')
            
        self._files.sort()
        self._steps = OrderedDict()
        for file in self._files:
            _, step, _ = file.split('.')
            step = int(step)
            if step in self._steps:
                self._steps[step].append(file)
            else:
                self._steps[step] = [file]
        # print("steps", self._steps.keys())

        # hacky way of figuring out gdims
        file = adios2py.file(os.path.join(path, next(iter(self._files))))
        assert len(file.vars) > 0
        var = next(iter(file.vars))
        self.gdims = np.asarray(file[var].shape)[0:3]
        file.close() # FIXME, should support with ... as 

        if L is not None:
            self.L = np.asarray(L)
            # FIXME, corner should also be passed (or rather read)
            self.corner = -.5 * self.L
        else:
            self.L = self.gdims
            self.corner = np.array([0., 0., 0.])

        self.x = np.linspace(self.corner[0], self.corner[0] + self.L[0], self.gdims[0], endpoint=False)
        self.y = np.linspace(self.corner[1], self.corner[1] + self.L[1], self.gdims[1], endpoint=False)
        self.z = np.linspace(self.corner[2], self.corner[2] + self.L[2], self.gdims[2], endpoint=False)
        
        first_step = next(iter(self._steps))
        self.activate_time(first_step)
        
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
        coords = { "x": self.x[start[0]:start[0]+count[0]],
                   "y": self.y[start[1]:start[1]+count[1]],
                   "z": self.z[start[2]:start[2]+count[2]], }
        return xr.DataArray(arr, dims=['x', 'y', 'z'], coords=coords)
                   

def read(run, time):
    py_e = run.read('py_he_e', [0, 0, 0], run.gdims)
    n_e = -run.read('rho_he_e', [0, 0, 0], run.gdims)
    return {"n_e" : n_e, "py_e" : py_e}