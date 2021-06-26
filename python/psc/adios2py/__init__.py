
import adios2
import numpy as np
import logging

_ad = adios2.ADIOS()

_dtype_map = {"float": np.float32, "double": np.float64}

class variable:
    def __init__(self, var, engine):
        self._var = var
        self._engine = engine
        self.name = var.Name()
        self.shape = var.Shape()[::-1]
        self.dtype = _dtype_map[var.Type()]
        
    def _set_selection(self, start, count):
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

        self._set_selection(sel_start, sel_count)

        arr = np.empty(arr_shape, dtype=self.dtype, order='F')
        print("reading ", self.name, sel_start, sel_count)
        self._engine.Get(self._var, arr, adios2.Mode.Sync)
        return arr

class file:
    def __init__(self, filename):
        logging.debug(f"adios2py: __init__ {filename}")
        self._io_name = f'io-{filename}'
        self._io = _ad.DeclareIO(self._io_name)
        self._engine = self._io.Open(filename, adios2.Mode.Read)
        self.vars = self._io.AvailableVariables().keys()
        
    def __enter__(self):
        logging.debug(f"adios2py: __enter__")
        return self
    
    def __exit__(self, type, value, traceback):
        logging.debug(f"adios2py: __exit__")
        self.close()
        
    def close(self):
        logging.debug(f"adios2py: close")
        self._engine.Close()
        _ad.RemoveIO(self._io_name)
        
    def __getitem__(self, varname):
        return variable(self._io.InquireVariable(varname), self._engine)

