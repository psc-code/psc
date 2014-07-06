
import numpy as np
from pylab import *

mx = 7
my = 7
mz = 7

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

def check_curr_3d(x, v, jx, jy, jz, rho_func):
    S0x = rho_func(x[0] - .5 * v[0])
    S1x = rho_func(x[0] + .5 * v[0])
    S0y = rho_func(x[1] - .5 * v[1])
    S1y = rho_func(x[1] + .5 * v[1])
    S0z = rho_func(x[2] - .5 * v[2])
    S1z = rho_func(x[2] + .5 * v[2])
    print "S0", S0x, S0y, S0z

    rho0, rho1 = zeros((mx, my, mz)), zeros((mx, my, mz))
    for k in xrange(mz):
        for j in xrange(my):
            for i in xrange(mx):
                rho0[i,j,k] = S0x[i] * S0y[j] * S0z[k]
                rho1[i,j,k] = S1x[i] * S1y[j] * S1z[k]

    drho = rho1 - rho0
    div_j = ((jx[1:  ,1:-1,1:-1] - jx[ :-1,1:-1,1:-1]) +
             (jy[1:-1,1:  ,1:-1] - jy[1:-1, :-1,1:-1]) + 
             (jz[1:-1,1:-1,1:  ] - jz[1:-1,1:-1, :-1]))
    err = drho[1:-1,1:-1,1:-1] + div_j
    if norm(err) > 1e-10:
        print "rho0"
        print rho0[1:-1,1:-1,1:-1].T
        print "rho1"
        print rho1[1:-1,1:-1,1:-1].T
        print "drho"
        print drho[1:-1,1:-1,1:-1].T
        print "div_j"
        print div_j.T
        print "err"
        print err.T
        raise Exception("continuity not satisfied")

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

def curr_2d_vb_split(jx, jy, i, xm, xp, dim=0):
    im = np.round(xm)
    i = i + im
    xm -= im
    xp -= im
    
    if dim == 2:
        dx = xp - xm
        print 'cell: i', i, 'xm', xm, 'xp', xp, 'dx', dx
        jx[i[0]  ,i[1]  ] += dx[0] * (.5 - xm[1] - .5 * dx[1])
        jx[i[0]  ,i[1]+1] += dx[0] * (.5 + xm[1] + .5 * dx[1])
        jy[i[0]  ,i[1]  ] += dx[1] * (.5 - xm[0] - .5 * dx[0])
        jy[i[0]+1,i[1]  ] += dx[1] * (.5 + xm[0] + .5 * dx[0])
        return
        
    if xp[dim] >= -.5 and xp[dim] <= .5:
        curr_2d_vb_split(jx, jy, i, xm, xp, dim + 1)
    else:
        idiff = zeros(2)
        idiff[dim] = sign(xp[dim])
        r = (.5 * idiff[dim] - xm[dim]) / (xp[dim] - xm[dim])
        xs = (1 - r) * xm + r * xp
        print 'splitting', xm, xs, xp
        curr_2d_vb_split(jx, jy, i, xm, xs, dim + 1)
        curr_2d_vb_split(jx, jy, i + idiff, xs - idiff, xp - idiff, dim + 1)
        

def curr_2d_vb(x, v, rho_func):
    x, v = array(x), array(v)
    jx = np.zeros((mx-1, my))
    jy = np.zeros((mx, my-1))

    xm = x - .5 * v
    xp = x + .5 * v

    xm -= .5
    xp -= .5

    im = zeros(2)
    curr_2d_vb_split(jx, jy, im, xm, xp)

    return jx, jy

def curr_3d_vb_split(jx, jy, jz, i, xm, xp, dim=0):
    im = np.round(xm)
    i = i + im
    xm -= im
    xp -= im
    
    if dim == 3:
        dx = xp - xm
        xa = .5 * (xm + xp)
        print 'cell: i', i, 'xm', xm, 'xp', xp, 'dx', dx
        h = dx[0] * dx[1] * dx[2] / 12.
        jx[i[0]  ,i[1]  ,i[2]  ] += dx[0] * (.5 - xa[1]) * (.5 - xa[2]) + h
        jx[i[0]  ,i[1]  ,i[2]+1] += dx[0] * (.5 - xa[1]) * (.5 + xa[2]) - h
        jx[i[0]  ,i[1]+1,i[2]  ] += dx[0] * (.5 + xa[1]) * (.5 - xa[2]) - h 
        jx[i[0]  ,i[1]+1,i[2]+1] += dx[0] * (.5 + xa[1]) * (.5 + xa[2]) + h
        jy[i[0]  ,i[1]  ,i[2]  ] += dx[1] * (.5 - xa[0]) * (.5 - xa[2]) + h
        jy[i[0]  ,i[1]  ,i[2]+1] += dx[1] * (.5 - xa[0]) * (.5 + xa[2]) - h
        jy[i[0]+1,i[1]  ,i[2]  ] += dx[1] * (.5 + xa[0]) * (.5 - xa[2]) - h
        jy[i[0]+1,i[1]  ,i[2]+1] += dx[1] * (.5 + xa[0]) * (.5 + xa[2]) + h
        jz[i[0]  ,i[1]  ,i[2]  ] += dx[2] * (.5 - xa[0]) * (.5 - xa[1]) + h
        jz[i[0]  ,i[1]+1,i[2]  ] += dx[2] * (.5 - xa[0]) * (.5 + xa[1]) - h
        jz[i[0]+1,i[1]  ,i[2]  ] += dx[2] * (.5 + xa[0]) * (.5 - xa[1]) - h
        jz[i[0]+1,i[1]+1,i[2]  ] += dx[2] * (.5 + xa[0]) * (.5 + xa[1]) + h
        return
        
    if xp[dim] >= -.5 and xp[dim] <= .5:
        curr_3d_vb_split(jx, jy, jz, i, xm, xp, dim + 1)
    else:
        idiff = zeros(3)
        idiff[dim] = sign(xp[dim])
        r = (.5 * idiff[dim] - xm[dim]) / (xp[dim] - xm[dim])
        xs = (1 - r) * xm + r * xp
        print 'splitting', xm, xs, xp
        curr_3d_vb_split(jx, jy, jz, i, xm, xs, dim + 1)
        curr_3d_vb_split(jx, jy, jz, i + idiff, xs - idiff, xp - idiff, dim + 1)
        

def curr_3d_vb(x, v, rho_func):
    x, v = array(x), array(v)
    jx = np.zeros((mx-1, my, mz))
    jy = np.zeros((mx, my-1, mz))
    jz = np.zeros((mx, my, mz-1))

    xm = x - .5 * v
    xp = x + .5 * v

    xm -= .5
    xp -= .5

    im = zeros(3)
    curr_3d_vb_split(jx, jy, jz, im, xm, xp)

    return jx, jy, jz

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

def test_3d(curr_3d, x, v):
    print 'x =', x, "v =", v
    jx, jy, jz = curr_3d(x, v, rho_1st)
    check_curr_3d(x, v, jx, jy, jz, rho_1st)
    for k in xrange(mz):
        for j in xrange(my):
            for i in xrange(mx-1):
                if abs(jx[i,j,k]) > 1e-15:
                    print '(%g,%d,%g) jx %g' % (i+.5, j, k, jx[i,j,k])
    for k in xrange(mz):
        for j in xrange(my-1):
            for i in xrange(mx):
                if abs(jy[i,j,k]) > 1e-15:
                    print '(%d,%g,%g) jy %g' % (i, j+.5, k, jy[i,j,k])
    for k in xrange(mz-1):
        for j in xrange(my):
            for i in xrange(mx):
                if abs(jz[i,j,k]) > 1e-15:
                    print '(%d,%g,%g) jz %g' % (i, j, k+.5, jz[i,j,k])
    # plot_rho(x, v)
    # plot_j(jx, jy)
    # savefig("vb1.pdf")

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
#test_2d(curr_2d_bv, (2.1, 2.3), (.3,-.9))
#test_2d(curr_2d_vb, (1.8, 2.0), (.8, .8))
test_3d(curr_3d_vb, (2.2, 2.4, 1.2), (.8, .8, .8))
#test_2d_all(curr_2d_bv)
