
import adios2
import numpy as np

_ad = adios2.ADIOS()

class variable:
    def __init__(self, var):
        self._var = var
        self.shape = var.Shape()[::-1]

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
        var = self._io.InquireVariable(varname)
        var.SetSelection((sel_start[::-1], sel_count[::-1]))
        arr = np.empty(sel_count, dtype='f', order='F')
        self._engine.Get(var, arr, adios2.Mode.Sync)
        return arr[:,:,:,0]
    
    def __getitem__(self, varname):
        return variable(self._io.InquireVariable(varname))

