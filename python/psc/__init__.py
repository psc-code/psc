
import numpy as np
import psc.adios2py
import adios2
import os
import xarray as xr

_ad = adios2.ADIOS()

class run:
    def __init__(self, path, L, gdims):
        # FIXME, corner should also be passed (or rather read)
        self.path = path
        self.L = np.asarray(L)
        self.gdims = np.asarray(gdims, dtype=int)
        self.dx = self.L / self.gdims
        self.corner = -.5 * self.L
        self.x = np.linspace(self.corner[0], self.corner[0] + self.L[0], self.gdims[0], endpoint=False)
        self.y = np.linspace(self.corner[1], self.corner[1] + self.L[1], self.gdims[1], endpoint=False)
        self.z = np.linspace(self.corner[2], self.corner[2] + self.L[2], self.gdims[2], endpoint=False)
        
class reader:
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

    def __init__(self, run, what, time):
        self._run = run
        self._what = what
        self._time = time
        self._file = adios2py.file(run.path, what, time)
        if what in ('pfd', 'tfd'):
            self._varname = 'jeh'
        elif what in ('pfd_moments', 'tfd_moments'):
            self._varname = "all_1st"
        else:
            raise(f'{what} not supported!')
        
    def read(self, fldname, start, count):
        m = self._to_index(fldname)
        sel_start = np.array([start[0], start[1], start[2], m])
        sel_count = np.array([count[0], count[1], count[2], 1])
        arr = self._file.read(self._varname, sel_start, sel_count)
        coords = { "x": self._run.x[start[0]:start[0]+count[0]],
                   "y": self._run.y[start[1]:start[1]+count[1]],
                   "z": self._run.z[start[2]:start[2]+count[2]], }
        return xr.DataArray(arr, dims=['x', 'y', 'z'], coords=coords)
        
    def _to_index(self, fldname):
        if self._varname == 'jeh':
            return self._jeh_to_index[fldname]
        elif self._varname == 'all_1st':
            return self._all_1st_to_index[fldname]
        else:
            assert False
            
def read(run, time):
    pfd = reader(run, "pfd_moments", time)
    py_e = pfd.read('py_he_e', [0, 0, 0], run.gdims)
    n_e = -pfd.read('rho_he_e', [0, 0, 0], run.gdims)
    return {"n_e" : n_e, "py_e" : py_e}