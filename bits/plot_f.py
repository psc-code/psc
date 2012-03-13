import time
import numpy as np
import matplotlib
matplotlib.use('TkAgg') # do this before importing pylab
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlim(0, 2*np.pi)
plt.ylim(-.6, .6)

def animate():
    line = None
    for ts in xrange(0, 2000, 10):
        filename = "prt.%06d_p000000.asc" % ts
        data = np.genfromtxt(filename, unpack=True)
        if line is None:
            line, = ax.plot(data[3,:],data[6,:],'.')
        else:
            line.set_xdata(data[3,:])
            line.set_ydata(data[6,:])
        fig.canvas.draw()

win = fig.canvas.manager.window
fig.canvas.manager.window.after(100, animate)
plt.show()

