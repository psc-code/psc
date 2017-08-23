
#include "psc_generic_c.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
do_genc_push_part_xz(int p, struct psc_fields *pf, particle_range_t prts)
{
#define S0X(off) s0x[off+2]
#define S0Z(off) s0z[off+2]
#define S1X(off) s1x[off+2]
#define S1Z(off) s1z[off+2]

  creal s0x[5] = {}, s0z[5] = {}, s1x[5], s1z[5];

  creal dt = ppsc->dt;
  creal xl = .5f * dt;
  creal zl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqxs = ppsc->patch[p].dx[0] * fnqs / dt;
  creal fnqzs = ppsc->patch[p].dx[2] * fnqs / dt;
  creal dxi = 1.f / ppsc->patch[p].dx[0];
  creal dzi = 1.f / ppsc->patch[p].dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);

    // x^n, p^n -> x^(n+.5), p^n

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->xi += vxi * xl;
    part->zi += vzi * zl;

    creal u = part->xi * dxi;
    creal w = part->zi * dzi;
    int lg1 = particle_real_nint(u);
    int lg3 = particle_real_nint(w);
    creal h1 = lg1-u;
    creal h3 = lg3-w;

    creal gmx=.5f*(.5f+h1)*(.5f+h1);
    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0x=.75f-h1*h1;
    creal g0z=.75f-h3*h3;
    creal g1x=.5f*(.5f-h1)*(.5f-h1);
    creal g1z=.5f*(.5f-h3)*(.5f-h3);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0X(-1) = gmx;
    S0X(+0) = g0x;
    S0X(+1) = g1x;
    S0Z(-1) = gmz;
    S0Z(+0) = g0z;
    S0Z(+1) = g1z;

    u = part->xi * dxi - .5f;
    w = part->zi * dzi - .5f;
    int lh1 = particle_real_nint(u);
    int lh3 = particle_real_nint(w);
    h1=lh1 - u;
    h3=lh3 - w;
    creal hmx=.5f*(.5f+h1)*(.5f+h1);
    creal hmz=.5f*(.5f+h3)*(.5f+h3);
    creal h0x=.75f-h1*h1;
    creal h0z=.75f-h3*h3;
    creal h1x=.5f*(.5f-h1)*(.5f-h1);
    creal h1z=.5f*(.5f-h3)*(.5f-h3);

    // FIELD INTERPOLATION

#define INTERPOLATE_FIELD(m, gx, gz)					\
    (gz##mz*(gx##mx*F3(pf, m, l##gx##1-1,0,l##gz##3-1) +			\
	     gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3-1) +			\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3-1)) +			\
     gz##0z*(gx##mx*F3(pf, m, l##gx##1-1,0,l##gz##3  ) +			\
	     gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3  ) +			\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3  )) +			\
     gz##1z*(gx##mx*F3(pf, m, l##gx##1-1,0,l##gz##3+1) +			\
	     gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3+1) +			\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3+1)))			\
      
    creal exq = INTERPOLATE_FIELD(EX, h, g);
    creal eyq = INTERPOLATE_FIELD(EY, g, g);
    creal ezq = INTERPOLATE_FIELD(EZ, g, h);

    creal hxq = INTERPOLATE_FIELD(HX, g, h);
    creal hyq = INTERPOLATE_FIELD(HY, h, h);
    creal hzq = INTERPOLATE_FIELD(HZ, h, g);

     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    creal dq = dqs * part->qni / part->mni;
    creal pxm = part->pxi + dq*exq;
    creal pym = part->pyi + dq*eyq;
    creal pzm = part->pzi + dq*ezq;

    root = dq / creal_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
    creal taux = hxq*root;
    creal tauy = hyq*root;
    creal tauz = hzq*root;

    creal tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    creal pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.f*taux*tauy+2.f*tauz)*pym + 
		(2.f*taux*tauz-2.f*tauy)*pzm)*tau;
    creal pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
		(1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.f*tauy*tauz+2.f*taux)*pzm)*tau;
    creal pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
		(2.f*tauy*tauz-2.f*taux)*pym +
		(1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
    
    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    vxi = part->pxi * root;
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->xi += vxi * xl;
    part->zi += vzi * zl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal xi = part->xi + vxi * xl;
    creal zi = part->zi + vzi * zl;

    u = xi * dxi;
    w = zi * dzi;
    int k1 = particle_real_nint(u);
    int k3 = particle_real_nint(w);
    h1 = k1 - u;
    h3 = k3 - w;

    for (int i = -2; i <= 2; i++) {
      S1X(i) = 0.f;
      S1Z(i) = 0.f;
    }

    S1X(k1-lg1-1) = .5f*(.5f+h1)*(.5f+h1);
    S1X(k1-lg1+0) = .75f-h1*h1;
    S1X(k1-lg1+1) = .5f*(.5f-h1)*(.5f-h1);
    S1Z(k3-lg3-1) = .5f*(.5f+h3)*(.5f+h3);
    S1Z(k3-lg3+0) = .75f-h3*h3;
    S1Z(k3-lg3+1) = .5f*(.5f-h3)*(.5f-h3);

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1X(i) -= S0X(i);
      S1Z(i) -= S0Z(i);
    }

    int l1min, l3min, l1max, l3max;
    
    if (k1 == lg1) {
      l1min = -1; l1max = +1;
    } else if (k1 == lg1 - 1) {
      l1min = -2; l1max = +1;
    } else { // (k1 == lg1 + 1)
      l1min = -1; l1max = +2;
    }

    if (k3 == lg3) {
      l3min = -1; l3max = +1;
    } else if (k3 == lg3 - 1) {
      l3min = -2; l3max = +1;
    } else { // (k3 == lg3 + 1)
      l3min = -1; l3max = +2;
    }

    creal fnqx = part->qni * part->wni * fnqxs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      creal jxh = 0.f;
      for (int l1 = l1min; l1 < l1max; l1++) {
	creal wx = S1X(l1) * (S0Z(l3) + .5f*S1Z(l3));
	jxh -= fnqx*wx;
	F3(pf, JXI, lg1+l1,0,lg3+l3) += jxh;
      }
    }

    creal fnqy = vyi * part->qni * part->wni * fnqs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l1 = l1min; l1 <= l1max; l1++) {
	creal wy = S0X(l1) * S0Z(l3)
	  + .5f * S1X(l1) * S0Z(l3)
	  + .5f * S0X(l1) * S1Z(l3)
	  + (1.f/3.f) * S1X(l1) * S1Z(l3);
	creal jyh = fnqy*wy;
	F3(pf, JYI, lg1+l1,0,lg3+l3) += jyh;
      }
    }

    creal fnqz = part->qni * part->wni * fnqzs;
    for (int l1 = l1min; l1 <= l1max; l1++) {
      creal jzh = 0.f;
      for (int l3 = l3min; l3 < l3max; l3++) {
	creal wz = S1Z(l3) * (S0X(l1) + .5f*S1X(l1));
	jzh -= fnqz*wz;
	F3(pf, JZI, lg1+l1,0,lg3+l3) += jzh;
      }
    }
  }
}

void
psc_push_particles_generic_c_push_mprts_xz(struct psc_push_particles *push,
					   struct psc_mparticles *mprts,
					   struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_xz", 1., 0, 0);
  }
  
  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    particle_range_t prts = particle_range_mprts(mprts, p);
    psc_fields_zero_range(flds, JXI, JXI + 3);
    do_genc_push_part_xz(p, flds, prts);
  }
  prof_stop(pr);
}

