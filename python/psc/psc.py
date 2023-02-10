
from . import adios2py

import numpy as np


class RunInfo:
    """Global information about the PSC run

    Currently stores domain info.
    TODO: Should also know about timestep, species, whatever...
    """

    def __init__(self, file, length=None, corner=None):
        assert len(file.variables) > 0
        var = next(iter(file.variables))
        self.gdims = np.asarray(file[var].shape)[0:3]

        maybe_length_attr = file._io.InquireAttribute("length")
        if maybe_length_attr:
            self.length = maybe_length_attr.Data()
            self.corner = file._io.InquireAttribute("corner").Data()
        elif length is not None:
            self.length = np.asarray(length)
            if corner is not None:
                self.corner = np.asarray(corner)
            else:
                self.corner = -.5 * self.length
        else:
            self.length = self.gdims
            self.corner = np.array([0., 0., 0.])

        self.x = np.linspace(
            self.corner[0], self.corner[0] + self.length[0], self.gdims[0], endpoint=False)
        self.y = np.linspace(
            self.corner[1], self.corner[1] + self.length[1], self.gdims[1], endpoint=False)
        self.z = np.linspace(
            self.corner[2], self.corner[2] + self.length[2], self.gdims[2], endpoint=False)

    def __repr__(self):
        return f"Psc(gdims={self.gdims}, length={self.length}, corner={self.corner})"


def FieldToComponent(species):
    map = dict()
    map['jeh'] = {'jx_ec': 0, 'jy_ec': 1, 'jz_ec': 2,
                  'ex_ec': 3, 'ey_ec': 4, 'ez_ec': 5,
                  'hx_fc': 6, 'hy_fc': 7, 'hz_fc': 8}
    map['dive'] = {'dive': 0}
    map['rho'] = {'rho': 0}
    map['d_rho'] = {'d_rho': 0}
    map['div_j'] = {'div_j': 0}

    map['all_1st'] = {}
    moments = ['rho', 'jx', 'jy', 'jz', 'px', 'py',
               'pz', 'txx', 'tyy', 'tzz', 'txy', 'tyz', 'tzx']
    for i_sp, sp in enumerate(species):
        for i_mom, mom in enumerate(moments):
            map['all_1st'][f"{mom}_{sp}"] = i_mom + 13 * i_sp

    return map
