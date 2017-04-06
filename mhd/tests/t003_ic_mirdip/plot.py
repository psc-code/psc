#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import viscid
from viscid.plot import mpl

steps = range(0, 1)

flds = ['divB', 'rr', 'bx', 'by', 'bz', 'vx', 'vy', 'vz', 'pp']

plot_kwargs = dict()

f = viscid.load_file("run.3d.xdmf")

for step in steps:
    print("Plotting step {}".format(step))
    f.activate_time(step)
    for fld in flds:
        mpl.plt.figure()

        dat = f[fld]
        # if fld == "divB":
        #     dat *= f["ymask"]
                                   
        mpl.plot(dat, "y=0.f", **plot_kwargs)

        mpl.plt.savefig("%s-xz-%06d.png" % (fld, step), dpi=200)
        mpl.plt.close()
