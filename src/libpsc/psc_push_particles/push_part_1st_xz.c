
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
do_push_part_1st_xz(int p, fields_t *pf, particles_t *pp)
{
#define S0X(off) s0x[off+1]
#define S0Z(off) s0z[off+1]
#define S1X(off) s1x[off+1]
#define S1Z(off) s1z[off+1]

  creal s0x[4] = {}, s0z[4] = {}, s1x[4], s1z[4];

  creal dt = ppsc->dt;
  creal xl = .5f * dt;
  creal zl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqxs = ppsc->dx[0] * fnqs / dt;
  creal fnqzs = ppsc->dx[2] * fnqs / dt;
  creal dxi = 1.f / ppsc->dx[0];
  creal dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    // x^n, p^n -> x^(n+.5), p^n

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->xi += vxi * xl;
    part->zi += vzi * zl;

    creal u = (part->xi - patch->xb[0]) * dxi;
    creal w = (part->zi - patch->xb[2]) * dzi;
    int lg1 = fint(u);
    int lg3 = fint(w);
    creal h1 = u - lg1;
    creal h3 = w - lg3;

    creal g0x = 1.f - h1;
    creal g0z = 1.f - h3;
    creal g1x = h1;
    creal g1z = h3;

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0X(+0) = g0x;
    S0X(+1) = g1x;
    S0Z(+0) = g0z;
    S0Z(+1) = g1z;

    u = (part->xi - patch->xb[0]) * dxi - .5f;
    w = (part->zi - patch->xb[2]) * dzi - .5f;
    int lh1 = fint(u);
    int lh3 = fint(w);
    h1 = u - lh1;
    h3 = w - lh3;
    creal h0x = 1.f - h1;
    creal h0z = 1.f - h3;
    creal h1x = h1;
    creal h1z = h3;

    // FIELD INTERPOLATION

#define INTERPOLATE_FIELD(m, gx, gz)					\
    (gz##0z*(gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3  ) +			\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3  )) +			\
     gz##1z*(gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3+1) +			\
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

    u = (xi - patch->xb[0]) * dxi;
    w = (zi - patch->xb[2]) * dzi;
    int k1 = fint(u);
    int k3 = fint(w);
    h1 = u - k1;
    h3 = w - k3;

    for (int i = -1; i <= 2; i++) {
      S1X(i) = 0.f;
      S1Z(i) = 0.f;
    }

    S1X(k1-lg1+0) = 1.f - h1;
    S1X(k1-lg1+1) = h1;
    S1Z(k3-lg3+0) = 1.f - h3;
    S1Z(k3-lg3+1) = h3;

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = 0; i <= 1; i++) {
      S1X(i) -= S0X(i);
      S1Z(i) -= S0Z(i);
    }

    int l1min, l3min, l1max, l3max;
    
    if (k1 == lg1) {
      l1min = 0; l1max = +1;
    } else if (k1 == lg1 - 1) {
      l1min = -1; l1max = +1;
    } else { // (k1 == lg1 + 1)
      l1min = 0; l1max = +2;
    }

    if (k3 == lg3) {
      l3min = 0; l3max = +1;
    } else if (k3 == lg3 - 1) {
      l3min = -1; l3max = +1;
    } else { // (k3 == lg3 + 1)
      l3min = 0; l3max = +2;
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
psc_push_particles_1st_push_xz(struct psc_push_particles *push,
			       mparticles_base_t *particles_base,
			       mfields_base_t *flds_base)
{
  mparticles_t particles;
  psc_mparticles_get_from(&particles, particles_base);
  mfields_t *flds = psc_mfields_get_from(EX, EX + 6, flds_base);

  static int pr;
  if (!pr) {
    pr = prof_register("1st_part_xz", 1., 0, 0);
  }
  prof_start(pr);
  psc_mfields_zero(flds, JXI);
  psc_mfields_zero(flds, JYI);
  psc_mfields_zero(flds, JZI);

  psc_foreach_patch(ppsc, p) {
    do_push_part_1st_xz(p, &flds->f[p], &particles.p[p]);
  }
  prof_stop(pr);

  psc_mfields_put_to(flds, JXI, JXI + 3, flds_base);
  psc_mparticles_put_to(&particles, particles_base);
}

