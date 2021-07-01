
import adios2
import numpy as np
import logging

_ad = adios2.ADIOS()

_dtype_map = {"float": np.float32, "double": np.float64, "uint64_t": np.uint64,
              "int64_t": np.int64, "int32_t": np.int32, "uint32_t": np.uint32}


class variable:
    def __init__(self, var, engine):
        self._var = var
        self._engine = engine
        self.name = self._name()
        self.shape = self._shape()
        self.dtype = self._dtype()
        logging.debug(f"variable __init__ var {var} engine {engine}")

    def close(self):
        logging.debug("adios2py.variable close")
        self._var = None
        self._engine = None

    def _assert_not_closed(self):
        if not self._var:
            raise ValueError("adios2py: variable is closed")

    def _set_selection(self, start, count):
        self._assert_not_closed()

        self._var.SetSelection((start[::-1], count[::-1]))

    def _shape(self):
        self._assert_not_closed()

        return self._var.Shape()[::-1]

    def _name(self):
        self._assert_not_closed()

        return self._var.Name()

    def _dtype(self):
        self._assert_not_closed()

        return _dtype_map[self._var.Type()]

    def __getitem__(self, args):
        self._assert_not_closed()

        if not isinstance(args, tuple):
            args = (args, )

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
                if idx < 0:
                    idx += shape[d]
                sel_start[d] = idx
                sel_count[d] = 1
                continue

            raise RuntimeError(f"invalid args to __getitem__: {args}")

        # print("A sel_start", sel_start, "sel_count",
        #       sel_count, "arr_shape", arr_shape)

        for d in range(len(args), len(shape)):
            sel_start[d] = 0
            sel_count[d] = shape[d]
            arr_shape.append(sel_count[d])

        # print("B sel_start", sel_start, "sel_count",
        #       sel_count, "arr_shape", arr_shape)

        self._set_selection(sel_start, sel_count)

        arr = np.empty(arr_shape, dtype=self.dtype, order='F')
        print("reading ", self.name, args)
        self._engine.Get(self._var, arr, adios2.Mode.Sync)
        return arr

    def __repr__(self):
        return f"adios2py.variable(name={self.name}, shape={self.shape}, dtype={self.dtype}"


class file:
    def __init__(self, filename, mode='r'):
        logging.debug(f"adios2py: __init__ {filename}")
        assert mode == 'r'
        self._io_name = f'io-{filename}'
        self._io = _ad.DeclareIO(self._io_name)
        self._engine = self._io.Open(filename, adios2.Mode.Read)
        self._open_vars = {}

        self.variables = self._io.AvailableVariables().keys()

    def __enter__(self):
        logging.debug(f"adios2py: __enter__")
        return self

    def __exit__(self, type, value, traceback):
        logging.debug(f"adios2py: __exit__")
        self.close()

    def __del__(self):
        logging.debug(f"adios2py: __del__")
        if self._engine:
            self.close()

    def close(self):
        logging.debug(f"adios2py: close")
        logging.debug(f"open vars {self._open_vars}")
        for varname, var in self._open_vars.items():
            var.close()

        self._engine.Close()
        self._engine = None

        _ad.RemoveIO(self._io_name)
        self._io = None
        self._io_name = None

    def __getitem__(self, varname):
        var = variable(self._io.InquireVariable(varname), self._engine)
        self._open_vars[varname] = var
        return var
