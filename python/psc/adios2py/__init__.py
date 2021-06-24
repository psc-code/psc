
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
        start = np.zeros_like(shape)
        count = np.zeros_like(shape)
        for d in range(len(args)):
            arg = args[d]
            assert isinstance(arg, slice)
            assert arg.step is None or arg.step == 1
            first = 0 if arg.start is None else arg.start
            last = shape[d] if arg.stop is None else arg.stop
            start[d] = first
            count[d] = last - first

        #print(f"d {d} start {start} count {count}")
        self.set_selection(start, count)

        arr = np.empty(count, dtype='f', order='F')
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
        
    def read(self, varname, sel_start, sel_count):
        var = self[varname]
        return var[sel_start[0]:sel_start[0]+sel_count[0],
                   sel_start[1]:sel_start[1]+sel_count[1],
                   sel_start[2]:sel_start[2]+sel_count[2],
                   sel_start[3]:sel_start[3]+sel_count[3]]
    
    def __getitem__(self, varname):
        return variable(self._io.InquireVariable(varname), self._engine)

