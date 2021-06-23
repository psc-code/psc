
import adios2
import numpy as np

_ad = adios2.ADIOS()

class variable:
    def __init__(self, var, engine):
        self._var = var
        self._engine = engine
        self.shape = var.Shape()[::-1]
        
    def set_selection(self, start, count):
        self._var.SetSelection((start[::-1], count[::-1]))

    def __getitem__(self, arg):
        count = self._var.Count()[::-1]
        arr = np.empty(count, dtype='f', order='F')
        self._engine.Get(self._var, arr, adios2.Mode.Sync)
        return arr[arg]

class file:
    def __init__(self, filename):
        self._io_name = f'io-{filename}'
        self._io = _ad.DeclareIO(self._io_name)
        self._engine = self._io.Open(filename, adios2.Mode.Read)
        self.vars = self._io.AvailableVariables().keys()
        
    def close(self):
        self._engine.Close()
        _ad.RemoveIO(self._io_name)
        
    def read(self, varname, sel_start, sel_count):
        var = self[varname]
        var.set_selection(sel_start, sel_count)
        return var[:]
    
    def __getitem__(self, varname):
        return variable(self._io.InquireVariable(varname), self._engine)

