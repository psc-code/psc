
import numpy as np
from pylab import *

colors = ['#ffc0c0', '#c0ffc0', '#c0c0ff', '#ffffc0']
colorss = ['#c0ffff', '#ffc0ff', '#ffffc0']
colors2 = ['r', 'g', 'b', '#b0b040']
colorss2 = ['#306060', '#603060', '#606030']
fontsize = 8

def plot_square(x, val):
    xp = x + .5
    xm = x - .5
    xc = array([xm[0], xp[0], xp[0], xm[0]])
    yc = array([xm[1], xm[1], xp[1], xp[1]])
    fill(xc, yc, colors[val])
    text(x[0], x[1], val, horizontalalignment="center",
         verticalalignment="center", size=fontsize)

def plot_square3(x, val):
    xp = x + .5
    xm = x - .5
    xc = array([xm[0], xp[0], xp[0], xm[0]])
    yc = array([xm[1], xm[1], xp[1], xp[1]])
    fill(xc, yc, colorss[val+1])
    text(x[0], x[1], val, horizontalalignment="center",
         verticalalignment="center", size=fontsize)

def plot_square2(x, val, j):
    xp = x + .5
    xm = x - .5
    xc = array([xm[0], xp[0], xp[0], xm[0]])
    yc = array([xm[1], xm[1], xp[1], xp[1]])
    fill(xc, yc, colors[j])
    text(x[0], x[1], val, horizontalalignment="center",
         verticalalignment="center", size=fontsize)

def plot_square4(x, val, j):
    xp = x + .5
    xm = x - .5
    xc = array([xm[0], xp[0], xp[0], xm[0]])
    yc = array([xm[1], xm[1], xp[1], xp[1]])
    fill(xc, yc, colorss[j])
    text(x[0], x[1], val, horizontalalignment="center",
         verticalalignment="center", size=fontsize)

def p1():
    ax = axes(frameon=False)
    blocks = [[0,0,1,0,0,1,0,0],
              [1, 1, 0, 1, 2, 1, 0],
              [1, 2, 2, 3, 1, 2],
              [2, 3, 3, 3, 2, 3]]

    text(-1, -0, "cell indices", horizontalalignment="right",
         verticalalignment="center", size=fontsize)
    ii = 0
    ib = 0
    for b in blocks:
        for i, n in enumerate(b):
            text(ii + i, 1, ib, horizontalalignment="center",
                 verticalalignment="center", size=fontsize)
            x = array([ii + i,0])
            plot_square(x, n)
            ib += 1
        ii += len(b) + 1

    nr_blocks = len(blocks)
    counts = zeros((nr_blocks, 4), dtype="int")
    for bn, b in enumerate(blocks):
        for n in b:
            counts[bn][n] += 1

    for j in xrange(4):
        text(-3, -(j+1.5), "# %d's" % j, horizontalalignment="center",
             verticalalignment="center", size=fontsize, color=colors2[j])
        
        ii = 0
        for bn, b in enumerate(blocks):
            text(ii, -(j+1.5), counts[bn][j], horizontalalignment="center",
                 verticalalignment="center", size=fontsize, color=colors2[j])
            if bn < nr_blocks - 1:
                arrow(ii+.5, -(j+6), len(b) - 1, 0, facecolor="#c0c0c0",
                      edgecolor='#c0c0c0', head_width=.3, head_length=1)
            elif j < 3:
                arrow(ii-.5, -(j+6), -ii+2, -1, facecolor="#c0c0c0",
                      edgecolor='#c0c0c0', head_width=.3, head_length=1)
            ii += len(b) + 1
            
    pfxsum = zeros((nr_blocks, 4), dtype="int")
    mysum = 0
    for j in xrange(4):
        for bn in xrange(nr_blocks):
            pfxsum[bn][j] = mysum
            mysum += counts[bn][j]

    text(-3, -6, "pfx sum", horizontalalignment="center",
         verticalalignment="center", size=fontsize)
        
    for j in xrange(4):
        ii = 0
        for bn, b in enumerate(blocks):
            text(ii, -(j+6), pfxsum[bn][j], horizontalalignment="center",
                 verticalalignment="center", size=fontsize, color=colors2[j])
            ii += len(b) + 1

    text(-1, -11, "scan cells", horizontalalignment="right",
         verticalalignment="center", size=fontsize)

    targets = blocks[:]
    new_blocks = zeros((27,),dtype="int")
    ii = 0
    for bn, b in enumerate(blocks):
        for i, n in enumerate(b):
            x = array([ii + i, -11])
            tgt = pfxsum[bn][n]
            targets[bn][i] = tgt
            new_blocks[tgt] = n
            plot_square2(x, tgt, n)
            pfxsum[bn][n] += 1
        ii += len(b) + 1

    text(-1, -12.5, "move", horizontalalignment="right",
         verticalalignment="center", size=fontsize)
    ii = 0
    ib = 0
    for bn, b in enumerate(blocks):
        for i, n in enumerate(b):
            x = array([ii + i, -12.5])
            nn = new_blocks[ib + i]
            plot_square(x, nn)
        ii += len(b) + 1
        ib += len(b)

    text(-1, -14, "cell bounds", horizontalalignment="right",
         verticalalignment="center", size=fontsize)
    offsets = array([0, 8, 16, 22, 27])
    ii = 0
    for b in xrange(4):
        for i in xrange(offsets[b], offsets[b+1]):
            x = array([ii, -14])
            nn = new_blocks[i]
            plot_square(x, nn)
            ii += 1
        ii += 1

    gca().set_aspect(True)
    gca().get_xaxis().set_visible(False)
    gca().get_yaxis().set_visible(False)
    savefig("sort.pdf")

def p2():
    ax = axes(frameon=False)
    blocks = [[0,0,1,0,0,1,0,0],
              [1, 1, 0, 1, 2, 1, 0],
              [1, 2, 2, 3, 1, 2],
              [2, 3, 3, 3, 2, 3]]
    offsets = [0, 8, 15, 21, 27]
    offsets2 = [0, 9, 17, 24, 31]

    text(-1, -0, "cell indices", horizontalalignment="right",
         verticalalignment="center", size=fontsize)
    ii = 0
    ib = 0
    for b in blocks:
        for i, n in enumerate(b):
            text(ii + i, 1, ib, horizontalalignment="center",
                 verticalalignment="center", size=fontsize)
            x = array([ii + i,0])
            plot_square(x, n)
            ib += 1
        ii += len(b) + 1

    for j, b in enumerate(blocks):
        for i, n in enumerate(b):
            b[i] = n - j

    text(-1, -1.5, "shifted", horizontalalignment="right",
         verticalalignment="center", size=fontsize)
    ii = 0
    for b in blocks:
        for i, n in enumerate(b):
            x = array([ii + i,-1.5])
            plot_square3(x, n)
        ii += len(b) + 1

    nr_blocks = len(blocks)
    counts = zeros((nr_blocks, 3), dtype="int")
    for bn, b in enumerate(blocks):
        for n in b:
            counts[bn][n+1] += 1

    for j in xrange(3):
        text(-3, -(j+3), "# %d's" % (j-1), horizontalalignment="center",
             verticalalignment="center", size=fontsize, color=colorss2[j])
        
        ii = 0
        for bn, b in enumerate(blocks):
            text(ii, -(j+3), counts[bn][j], horizontalalignment="center",
                 verticalalignment="center", size=fontsize, color=colorss2[j])
            ii += len(b) + 1
            
    pfxsum = zeros((nr_blocks, 3), dtype="int")
    mysum = 0
    last = None
    for bn in xrange(1, nr_blocks+1):
        for j in xrange(-2,1):
            if bn + j < 0 or bn + j >= nr_blocks:
                continue
            if last is not None:
                arrow(offsets2[last[0]], - last[1] - 7,
                      .85*(offsets2[(bn + j)] - offsets2[last[0]]), .85*(-(-j - last[1])),
                      facecolor="#c0c0c0",
                      edgecolor='#c0c0c0', head_width=.3, head_length=1)
            pfxsum[bn + j, -j] = mysum
            mysum += counts[bn + j, -j]
            last = (bn + j, -j)

    text(-3, -7, "pfx sum", horizontalalignment="center",
         verticalalignment="center", size=fontsize)
        
    for j in xrange(3):
        ii = 0
        for bn, b in enumerate(blocks):
            text(ii, -(j+7), pfxsum[bn][j], horizontalalignment="center",
                 verticalalignment="center", size=fontsize, color=colorss2[j])
            ii += len(b) + 1

    text(-1, -11, "scan cells", horizontalalignment="right",
         verticalalignment="center", size=fontsize)

    targets = blocks[:]
    new_blocks = zeros((27,),dtype="int")
    ii = 0
    for bn, b in enumerate(blocks):
        for i, n in enumerate(b):
            x = array([ii + i, -11])
            tgt = pfxsum[bn][n+1]
            targets[bn][i] = tgt
            new_blocks[tgt] = n + bn
            plot_square4(x, tgt, n+1)
            pfxsum[bn][n+1] += 1
        ii += len(b) + 1

    text(-1, -12.5, "move and shift back", horizontalalignment="right",
         verticalalignment="center", size=fontsize)
    ii = 0
    ib = 0
    for bn, b in enumerate(blocks):
        for i, n in enumerate(b):
            x = array([ii + i, -12.5])
            nn = new_blocks[ib + i]
            print nn
            plot_square(x, nn)
        ii += len(b) + 1
        ib += len(b)

    text(-1, -14, "cell bounds", horizontalalignment="right",
         verticalalignment="center", size=fontsize)
    offsets = array([0, 8, 16, 22, 27])
    ii = 0
    for b in xrange(4):
        for i in xrange(offsets[b], offsets[b+1]):
            x = array([ii, -14])
            nn = new_blocks[i]
            plot_square(x, nn)
            ii += 1
        ii += 1

    gca().set_aspect(True)
    gca().get_xaxis().set_visible(False)
    gca().get_yaxis().set_visible(False)
    savefig("sort2.pdf")

figure()
p1()
figure()
p2()
print "done"
