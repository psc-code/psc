
#include "psc_generic_c.h"

#include "util/profile.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void
genc_push_part_xz()
{
  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_xz", 1., 0, psc.n_part * 12 * sizeof(real));
  }
  prof_start(pr);
 
#define S0X(off) s0x[off+1]
#define S0Z(off) s0z[off+1]
#define S1X(off) s1x[off+2]
#define S1Z(off) s1z[off+2]

  real s0x[3], s0z[3], s1x[5], s1z[5];

  struct psc_genc *genc = psc.c_ctx;

  real dt = psc.dt;
  real xl = .5f * dt;
  real zl = .5f * dt;
  real dqs = .5f * psc.coeff.eta * dt;
  real fnqs = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
  real fnqxs = psc.dx[0] * fnqs / dt;
  real fnqzs = psc.dx[2] * fnqs / dt;
  real dxi = 1. / psc.dx[0];
  real dyi = 1. / psc.dx[1];
  real dzi = 1. / psc.dx[2];

  memset(&F3(JXI, psc.ilg[0], psc.ilg[1], psc.ilg[2]), 0,
	 psc.fld_size * sizeof(f_real));
  memset(&F3(JYI, psc.ilg[0], psc.ilg[1], psc.ilg[2]), 0,
	 psc.fld_size * sizeof(f_real));
  memset(&F3(JZI, psc.ilg[0], psc.ilg[1], psc.ilg[2]), 0,
	 psc.fld_size * sizeof(f_real));
  
  for (int n = 0; n < psc.n_part; n++) {
    struct c_particle *part = &genc->part[n];

    // x^n, p^n -> x^(n+0.5), p^n

    real root = 1. / sqrt(1. + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    real vxi = part->pxi * root;
    real vyi = part->pxi * root;
    real vzi = part->pzi * root;

    part->xi += vxi * xl;
    part->zi += vzi * zl;

    real u = part->xi * dxi;
    real v = part->yi * dyi;
    real w = part->zi * dzi;
    int j1 = nint(u);
    int j2 = nint(v);
    int j3 = nint(w);
    real h1 = j1-u;
    real h2 = j2-v;
    real h3 = j3-w;

    real gmx=0.5*(0.5+h1)*(0.5+h1);
    real gmz=0.5*(0.5+h3)*(0.5+h3);
    real g0x=0.75-h1*h1;
    real g0z=0.75-h3*h3;
    real g1x=0.5*(0.5-h1)*(0.5-h1);
    real g1z=0.5*(0.5-h3)*(0.5-h3);

    // CHARGE DENSITY FORM FACTOR AT (n+0.5)*dt 

    S0X(-1) = 0.5*(1.5-fabs(h1-1.0))*(1.5-fabs(h1-1.0));
    S0X(+0) = 0.75-fabs(h1)*fabs(h1);
    S0X(+1) = 0.5*(1.5-fabs(h1+1.0))*(1.5-fabs(h1+1.0));
    S0Z(-1) = 0.5*(1.5-fabs(h3-1.0))*(1.5-fabs(h3-1.0));
    S0Z(+0) = 0.75-fabs(h3)*fabs(h3);
    S0Z(+1) = 0.5*(1.5-fabs(h3+1.0))*(1.5-fabs(h3+1.0));

    u = part->xi*dxi-0.5;
    v = part->yi*dyi;
    w = part->zi*dzi-0.5;
    int l1=nint(u);
    int l2=nint(v);
    int l3=nint(w);
    h1=l1-u;
    h2=l2-v;
    h3=l3-w;
    real hmx=0.5*(0.5+h1)*(0.5+h1);
    real hmz=0.5*(0.5+h3)*(0.5+h3);
    real h0x=0.75-h1*h1;
    real h0z=0.75-h3*h3;
    real h1x=0.5*(0.5-h1)*(0.5-h1);
    real h1z=0.5*(0.5-h3)*(0.5-h3);

    // FIELD INTERPOLATION

    real exq=gmz*(hmx*F3(EX, l1-1,j2,j3-1)
		  +h0x*F3(EX, l1,j2,j3-1)
		  +h1x*F3(EX, l1+1,j2,j3-1))
      +g0z*(hmx*F3(EX, l1-1,j2,j3)
	    +h0x*F3(EX, l1,j2,j3)
	    +h1x*F3(EX, l1+1,j2,j3))
      +g1z*(hmx*F3(EX, l1-1,j2,j3+1)
	    +h0x*F3(EX, l1,j2,j3+1)
	    +h1x*F3(EX, l1+1,j2,j3+1));

    real eyq=gmz*(gmx*F3(EY, j1-1,l2,j3-1)
                   +g0x*F3(EY, j1,l2,j3-1)
                   +g1x*F3(EY, j1+1,l2,j3-1))
              +g0z*(gmx*F3(EY, j1-1,l2,j3)
                   +g0x*F3(EY, j1,l2,j3)
                   +g1x*F3(EY, j1+1,l2,j3))
              +g1z*(gmx*F3(EY, j1-1,l2,j3+1)
                   +g0x*F3(EY, j1,l2,j3+1)
		    +g1x*F3(EY, j1+1,l2,j3+1));

    real ezq=hmz*(gmx*F3(EZ, j1-1,j2,l3-1)
		   +g0x*F3(EZ, j1,j2,l3-1)
                   +g1x*F3(EZ, j1+1,j2,l3-1))
              +h0z*(gmx*F3(EZ, j1-1,j2,l3)
                   +g0x*F3(EZ, j1,j2,l3)
                   +g1x*F3(EZ, j1+1,j2,l3))
              +h1z*(gmx*F3(EZ, j1-1,j2,l3+1)
                   +g0x*F3(EZ, j1,j2,l3+1)
		    +g1x*F3(EZ, j1+1,j2,l3+1));

    real bxq=hmz*(gmx*F3(BX, j1-1,l2,l3-1)
                   +g0x*F3(BX, j1,l2,l3-1)
                   +g1x*F3(BX, j1+1,l2,l3-1))
              +h0z*(gmx*F3(BX, j1-1,l2,l3)
                   +g0x*F3(BX, j1,l2,l3)
                   +g1x*F3(BX, j1+1,l2,l3))
              +h1z*(gmx*F3(BX, j1-1,l2,l3+1)
                   +g0x*F3(BX, j1,l2,l3+1)
		    +g1x*F3(BX, j1+1,l2,l3+1));

    real byq=hmz*(hmx*F3(BY, l1-1,j2,l3-1)
                   +h0x*F3(BY, l1,j2,l3-1)
                   +h1x*F3(BY, l1+1,j2,l3-1))
              +h0z*(hmx*F3(BY, l1-1,j2,l3)
                   +h0x*F3(BY, l1,j2,l3)
                   +h1x*F3(BY, l1+1,j2,l3))
              +h1z*(hmx*F3(BY, l1-1,j2,l3+1)
                   +h0x*F3(BY, l1,j2,l3+1)
		    +h1x*F3(BY, l1+1,j2,l3+1));

    real bzq=gmz*(hmx*F3(BZ, l1-1,l2,j3-1)
		    +h0x*F3(BZ, l1,l2,j3-1)
                   +h1x*F3(BZ, l1+1,l2,j3-1))
              +g0z*(hmx*F3(BZ, l1-1,l2,j3)
                   +h0x*F3(BZ, l1,l2,j3)
                   +h1x*F3(BZ, l1+1,l2,j3))
              +g1z*(hmx*F3(BZ, l1-1,l2,j3+1)
                   +h0x*F3(BZ, l1,l2,j3+1)
		    +h1x*F3(BZ, l1+1,l2,j3+1));

     // c x^(n+0.5), p^n -> x^(n+1.0), p^(n+1.0) 

    real dq = dqs * part->qni / part->mni;
    real pxm = part->pxi + dq*exq;
    real pym = part->pyi + dq*eyq;
    real pzm = part->pzi + dq*ezq;

    root = dq / sqrt(1. + pxm*pxm + pym*pym + pzm*pzm);
    real taux = bxq*root;
    real tauy = byq*root;
    real tauz = bzq*root;

    real tau = 1. / (1. + taux*taux + tauy*tauy + tauz*tauz);
    real pxp = ((1.0+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.0*taux*tauy+2.0*tauz)*pym + 
		(2.0*taux*tauz-2.0*tauy)*pzm)*tau;
    real pyp = ((2.0*taux*tauy-2.0*tauz)*pxm +
		(1.0-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.0*tauy*tauz+2.0*taux)*pzm)*tau;
    real pzp = ((2.0*taux*tauz+2.0*tauy)*pxm +
		(2.0*tauy*tauz-2.0*taux)*pym +
		(1.0-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
    
    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    root = 1. / sqrt(1. + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    vxi = part->pxi * root;
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->xi += vxi * xl;
    part->zi += vzi * zl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)

    real xi = part->xi + vxi * xl;
    real zi = part->zi + vzi * zl;

    u = xi * dxi;
    w = zi * dzi;
    int k1 = nint(u);
    int k3 = nint(w);
    h1 = k1 - u;
    h3 = k3 - w;

    for (int i = -2; i <= 2; i++) {
      S1X(i) = 0.;
      S1Z(i) = 0.;
    }

    S1X(k1-j1-1) = 0.5*(1.5-fabs(h1-1.0))*(1.5-fabs(h1-1.0));
    S1X(k1-j1+0) = 0.75-fabs(h1)*fabs(h1);
    S1X(k1-j1+1) = 0.5*(1.5-fabs(h1+1.0))*(1.5-fabs(h1+1.0));
    S1Z(k3-j3-1) = 0.5*(1.5-fabs(h3-1.0))*(1.5-fabs(h3-1.0));
    S1Z(k3-j3+0) = 0.75-fabs(h3)*fabs(h3);
    S1Z(k3-j3+1) = 0.5*(1.5-fabs(h3+1.0))*(1.5-fabs(h3+1.0));

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1X(i) -= S0X(i);
      S1Z(i) -= S0Z(i);
    }

    int l1min, l3min, l1max, l3max;
    
    if (k1 == j1) {
      l1min = -1; l1max = +1;
    } else if (k1 == j1 - 1) {
      l1min = -2; l1max = +1;
    } else { // (k1 == j1 + 1)
      l1min = -1; l1max = +2;
    }

    if (k3 == j3) {
      l3min = -1; l3max = +1;
    } else if (k3 == j3 - 1) {
      l3min = -2; l3max = +1;
    } else { // (k3 == j3 + 1)
      l3min = -1; l3max = +2;
    }

    real jxh[6][5][5] = {};
    real jyh[5][6][5] = {};
    real jzh[5][5][6] = {};

#define JXH(i,j,k) jxh[i+3][j+2][k+2]
#define JYH(i,j,k) jyh[i+2][j+3][k+2]
#define JZH(i,j,k) jzh[i+2][j+2][k+3]

    real fnqx = part->qni * part->wni * fnqxs;
    real fnqy = vyi * part->qni * part->wni * fnqs;
    real fnqz = part->qni * part->wni * fnqzs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l1 = l1min; l1 <= l1max; l1++) {
	real wx = S1X(l1) * (S0Z(l3) + .5*S1Z(l3));
	real wy = S0X(l1) * S0Z(l3)
	  + .5 * S1X(l1) * S0Z(l3)
	  + .5 * S0X(l1) * S1Z(l3)
	  + (1./3.) * S1X(l1) * S1Z(l3);
	real wz = S1Z(l3) * (S0X(l1) + .5*S1X(l1));

	JXH(l1,0,l3) = JXH(l1-1,0,l3)-fnqx*wx;
	JYH(l1,0,l3) = fnqy*wy;
	JZH(l1,0,l3) = JZH(l1,0,l3-1)-fnqz*wz;

	F3(JXI, j1+l1,j2,j3+l3) += JXH(l1,0,l3);
	F3(JYI, j1+l1,j2,j3+l3) += JYH(l1,0,l3);
	F3(JZI, j1+l1,j2,j3+l3) += JZH(l1,0,l3);
      }
    }

  }

  prof_stop(pr);
}
