
import numpy as np
from pylab import *

mx = 7
my = 7

def nint(x):
    return int(x + .5 + 10) - 10

def fint(x):
    return int(x + 10) - 10

def rho_1st(x):
    ix = fint(x)
    h = x - ix
    S = np.zeros(mx)
    S[ix  ] = 1 - h
    S[ix+1] = h
    return S

def rho_2nd(x):
    ix = nint(x)
    h = x - ix
    S = np.zeros(mx)
    S[ix-1] = .5 * (.5 - h)**2
    S[ix  ] = .75 - h**2
    S[ix+1] = .5 * (.5 + h)**2
    return S

def check_curr_2d(x, v, jx, jy, rho_func):
    S0x = rho_func(x[0] - .5 * v[0])
    S1x = rho_func(x[0] + .5 * v[0])
    S0y = rho_func(x[1] - .5 * v[1])
    S1y = rho_func(x[1] + .5 * v[1])

    rho0, rho1 = zeros((mx, my)), zeros((mx, my))
    for j in xrange(my):
        for i in xrange(mx):
            rho0[i,j] = S0x[i] * S0y[j]
            rho1[i,j] = S1x[i] * S1y[j]

    drho = rho1 - rho0
    div_j = (jx[1:,1:-1] - jx[:-1,1:-1]) + (jy[1:-1,1:] - jy[1:-1,:-1])
    err = drho[1:-1,1:-1] + div_j
    if norm(err) > 1e-10:
        print "drho"
        print drho[1:-1,1:-1].T
        print "div_j"
        print div_j.T
        print "err"
        print err.T
        raise Exception("continuity not satisfied")
    for  j in xrange(my-1):
        for i in xrange(mx-1):
            pass
    

def curr_2d_psc(x, v, rho_func):
    S0x = rho_func(x[0] - .5 * v[0])
    S1x = rho_func(x[0] + .5 * v[0])
    S0y = rho_func(x[1] - .5 * v[1])
    S1y = rho_func(x[1] + .5 * v[1])

    jx = np.zeros((mx-1, my))
    for j in xrange(my):
        jxh = 0
        for i in xrange(mx-1):
            jxh -= (S1x[i] - S0x[i]) * .5 * (S0y[j] + S1y[j])
            jx[i,j] = jxh

    jy = np.zeros((mx, my-1))
    for i in xrange(mx):
        jyh = 0
        for j in xrange(my-1):
            jyh -= (S1y[j] - S0y[j]) * .5 * (S0x[i] + S1x[i])
            jy[i,j] = jyh

    return jx, jy

def curr_2d_bv_cell(jx, jy, i, x, dx, dxt=None, off=None):
    jx[i[0]  ,i[1]  ] += dx[0] * (.5 - x[1] - .5 * dx[1])
    jx[i[0]  ,i[1]+1] += dx[0] * (.5 + x[1] + .5 * dx[1])
    jy[i[0]  ,i[1]  ] += dx[1] * (.5 - x[0] - .5 * dx[0])
    jy[i[0]+1,i[1]  ] += dx[1] * (.5 + x[0] + .5 * dx[0])
    # print "i %d %d jx %g" % (i[0], i[1], jx[i[0], i[1]])
    # print "i %d %d jx %g" % (i[0], i[1]+1, jx[i[0], i[1]+1])
    # print "i %d %d jy %g" % (i[0], i[1], jy[i[0], i[1]])
    # print "i %d %d jy %g" % (i[0], i[1], jy[i[0]+1, i[1]])
    if dxt is not None:
        dxt -= dx
        x += dx
    if off is not None:
        i += off
        x -= off

def calc_dx1(dx1, x, dx, off):
    if off[1] == 0:
        dx1[0] = .5 * off[0] - x[0]
        dx1[1] = dx[1] / dx[0] * dx1[0]
    else:
        dx1[1] = .5 * off[1] - x[1]
        dx1[0] = dx[0] / dx[1] * dx1[1]
        
def curr_2d_bv(x, v, rho_func):
    x, v= array(x), array(v)
    jx = np.zeros((mx-1, my))
    jy = np.zeros((mx, my-1))

    xm = x - .5 * v
    xp = x + .5 * v
    im = array((fint(xm[0]), fint(xm[1])))
    ip = array((fint(xp[0]), fint(xp[1])))

    i = im
    idiff = ip - im
    dx = xp - xm
    x = xm - (i + .5)
    
    dx1 = empty(2)
    off = empty(2)
    second_dir = -1
    # FIXME, make sure we never div-by-zero?
    if idiff[0] == 0 and idiff[1] == 0:
        first_dir = -1
    elif idiff[0] == 0:
        first_dir = 1
    elif idiff[1] == 0:
        first_dir = 0
    else:
        dx1[0] = .5 * idiff[0] - x[0]
        dx1[1] = dx[1] / dx[0] * dx1[0]
        if abs(x[1] + dx1[1]) > .5:
            first_dir = 1
        else:
            first_dir = 0
        second_dir = 1 - first_dir

    if first_dir >= 0:
        off[1-first_dir] = 0
        off[first_dir] = idiff[first_dir]
        calc_dx1(dx1, x, dx, off)
        curr_2d_bv_cell(jx, jy, i, x, dx1, dx, off)

    if second_dir >= 0:
        off[first_dir] = 0
        off[second_dir] = idiff[second_dir]
        calc_dx1(dx1, x, dx, off)
        curr_2d_bv_cell(jx, jy, i, x, dx1, dx, off)

    curr_2d_bv_cell(jx, jy, i, x, dx)
        
    return jx, jy

def plot_cloud(x, **kwargs):
    xp = x + .5
    xm = x - .5
    plot([xm[0],xp[0],xp[0],xm[0],xm[0]], [xm[1], xm[1], xp[1], xp[1], xm[1]], *kwargs)
    i, j = int(x[0]), int(x[1])
    plot([i,i,i+1,i+1],[j,j+1,j,j+1],'o')

def plot_rho(x, v):
    x = array(x)
    v = array(v)
    plot_cloud(x - .5 * v)
    plot_cloud(x + .5 * v)
    xlim(0, 4)
    ylim(0, 4)
    grid(True)
    gca().set_aspect(True)

def plot_j(jx, jy):
    for i in xrange(mx-1):
        for j in xrange(my):
            if abs(jx[i,j]) > 1e-15:
                plot([i+.5, i+.5], [j-.5, j+.5], 'g')
    for j in xrange(my-1):
        for i in xrange(mx):
            if abs(jy[i,j]) > 1e-15:
                plot([i-.5, i+.5], [j+.5, j+.5], 'g')

def test_1d():
    for x in [2., 2.05, 2.45, 2.5, 2.55, 2.95]:
        print 'x =', x
        jx = curr_1d(x, 1, rho_1st)
        for i in xrange(mx-1):
            if abs(jx[i]) > 1e-15:
                print i+.5, 'jx', jx[i]
        print

def test_2d(curr_2d, x, v):
    #    for xv in [(2.0, 2.0), (2.05, 2.05), (2.45, 2.45), (2.5, 2.5), (2.95, 2.95)]:
    # (2,2.1) -> (3  , 2.5)
    # (2.25, 2.2) : (2  , 2.1) -> (2.5, 2.3)
    # (2.5,2) jx 0.4
    # (2.5,3) jx 0.1
    # (2,2.5) jy 0.15
    # (3,2.5) jy 0.05
    # (2.75, 2.4) : (2.5, 2.3) -> (3  , 2.5)
    # (2.5,2) jx 0.3
    # (2.5,3) jx 0.2
    # (2,2.5) jy 0.05
    # (3,2.5) jy 0.15
    # (2.25, 2.2) -> (3.25, 2.6):
    #(2.5,2) jx 0.4875
    #(3.5,2) jx 0.1125
    #(2.5,3) jx 0.2625
    #(3.5,3) jx 0.1375
    #(2,2.5) jy 0.1125
    #(3,2.5) jy 0.275
    #(4,2.5) jy 0.0125

    print 'x =', x, "v =", v
    jx, jy = curr_2d(x, v, rho_1st)
    check_curr_2d(x, v, jx, jy, rho_1st)
    for j in xrange(my):
        for i in xrange(mx-1):
            if abs(jx[i,j]) > 1e-15:
                print '(%g,%d) jx %g' % (i+.5, j, jx[i,j])
    for j in xrange(my-1):
        for i in xrange(mx):
            if abs(jy[i,j]) > 1e-15:
                print '(%d,%g) jy %g' % (i, j+.5, jy[i,j])
    plot_rho(x, v)
    plot_j(jx, jy)
    savefig("vb1.pdf")
    show()

def test_2d_all(curr_2d):
    fig = 0

    x = empty(2)
    v = empty(2)

    for v[1] in arange(-1., 1., 0.2):
        for v[0] in arange(-1., 1., 0.2):
            for x[1] in arange(2., 3., 0.1):
                for x[0] in arange(2., 3., 0.1):
                    try:
                        jx, jy = curr_2d(x, v, rho_1st)
                        check_curr_2d(x, v, jx, jy, rho_1st)
                    except:
                        print "x =", x, "v =", v
                        raise

# test_1d()
test_2d(curr_2d_bv, (2.1, 2.3), (.3,-.9))
#test_2d_all(curr_2d_bv)
