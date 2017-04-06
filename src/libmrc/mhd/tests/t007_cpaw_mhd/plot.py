#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import viscid
from viscid.plot import mpl

steps = range(0, 21)

flds = ['divB', 'rr', 'bx', 'by', 'bz', 'vx', 'vy', 'vz', 'pp']
flds = ['vy']

f = viscid.load_file("run.3d.xdmf")

for step in steps:
    print("Plotting step {}".format(step))
    f.activate_time(step)
    for fld in flds:
        mpl.plt.figure()

        plot_kwargs = dict()
        dat = f[fld]
        # if fld == "divB":
        #     dat *= f["ymask"]

        if fld == "vx":
            plot_kwargs["cmap"] = "hot"
            plot_kwargs["clim"] = (-150, 150)
                                   
        mpl.plot(dat, "y=0.f", **plot_kwargs)
        mpl.plt.ylim(-1.5e-3, 1.5e-3)

        mpl.plt.savefig("%s-xz-%06d.png" % (fld, step), dpi=200)
        mpl.plt.close()
