#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import viscid
from viscid.plot import mpl

steps = range(0, 101, 1)

flds = ['rr', 'bz']

plot_kwargs = dict(cmap='jet')

f = viscid.load_file("run.3d.xdmf")

for step in steps:
    print("Plotting step {}".format(step))
    f.activate_time(step)
    for fld in flds:
        mpl.plt.figure()

        dat = f[fld]
        if fld == "divB":
            dat *= f["ymask"]
                                   
        mpl.plot(dat, "y=0.f", **plot_kwargs)

        mpl.plt.savefig("%s-xz-%06d.png" % (fld, step), dpi=100)
        mpl.plt.close()
