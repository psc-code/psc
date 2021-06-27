"""
Eventually, we want to have a way to deal with an entire run at once, though right
now this is pretty much just out-of-date hacky stuff.
"""

import os

from .psc import RunInfo, FieldToComponent
from . import adios2py

import xarray as xr

from collections import namedtuple

field_entry = namedtuple('field_entry', ['filename', 'varname', 'index'])

class Run:
    def __init__(self, path, L=None, pfx='pfd'):
        self._path = path
        self._pfx = pfx

        # open first var in first file to figure out global dims
        self._files = [f for f in os.listdir(path) if f.startswith(self._pfx) and f.endswith('.bp')]
        if not self._files:
            raise RuntimeError(f"No {self._pfx} data files found in {path}")
            
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

        self.psc = RunInfo(os.path.join(path, next(iter(self._files))), length=L)

        self._fields_to_index = FieldToComponent(['he_e', 'e', 'i'])

        first_step = next(iter(self._steps))
        self.activate_time(first_step)
        
    def activate_time(self, step):
        if not step in self._steps:
            raise ValueError(f'step {step} not found!')
        self._step = step
        filenames = self._steps[step]

        self.fields = {}
        for filename in filenames:
            with adios2py.file(os.path.join(self._path, filename)) as file:
                for varname in file.variables:
                    fields_to_index = self._fields_to_index[varname]
                    # assert # of comps (and gdims?) match
                    for f in fields_to_index:
                        self.fields[f] = field_entry(filename=filename, varname=varname, index=fields_to_index[f])
        # print(f'fields {self.fields.keys()}')
        
    def read(self, fldname, start, count):
        field = self.fields[fldname]
        with adios2py.file(os.path.join(self._path, field.filename)) as file:
            var = file[field.varname]
            m = field.index
            arr = var[start[0]:start[0]+count[0],
                      start[1]:start[1]+count[1],
                      start[2]:start[2]+count[2],
                      m]

        coords = { "x": self.psc.x[start[0]:start[0]+count[0]],
                   "y": self.psc.y[start[1]:start[1]+count[1]],
                   "z": self.psc.z[start[2]:start[2]+count[2]], }
        return xr.DataArray(arr, dims=['x', 'y', 'z'], coords=coords)

