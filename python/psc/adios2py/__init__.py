
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

    def __getitem__(self, args):
        #print('args', args)
        shape = self.shape
        sel_start = np.zeros_like(shape)
        sel_count = np.zeros_like(shape)
        arr_shape = []
        for d, arg in enumerate(args):
            if isinstance(arg, slice):
                start, stop, step = arg.indices(shape[d])
                assert stop > start
                assert step == 1
                sel_start[d] = start
                sel_count[d] = stop - start
                arr_shape.append(sel_count[d])
                continue
                
            try:
                idx = int(arg)
            except:
                pass
            else:
                if idx < 0: idx += shape[d]
                sel_start[d] = idx
                sel_count[d] = 1
                continue
            
            raise RuntimeError(f"invalid args to __getitem__: {args}")

        self.set_selection(sel_start, sel_count)

        arr = np.empty(arr_shape, dtype='f', order='F')
        self._engine.Get(self._var, arr, adios2.Mode.Sync)
        return arr

class file:
    def __init__(self, filename):
        self._io_name = f'io-{filename}'
        self._io = _ad.DeclareIO(self._io_name)
        self._engine = self._io.Open(filename, adios2.Mode.Read)
        self.vars = self._io.AvailableVariables().keys()
        
    def close(self):
        self._engine.Close()
        _ad.RemoveIO(self._io_name)
        
    def __getitem__(self, varname):
        return variable(self._io.InquireVariable(varname), self._engine)
