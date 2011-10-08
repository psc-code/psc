
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
do_push_part_1st_yz(int p, fields_t *pf, particles_t *pp)
{
#define S0Y(off) s0y[off+1]
#define S0Z(off) s0z[off+1]
#define S1Y(off) s1y[off+1]
#define S1Z(off) s1z[off+1]

  creal s0y[4] = {}, s0z[4] = {}, s1y[4], s1z[4];

  creal dt = ppsc->dt;
  creal yl = .5f * dt;
  creal zl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqys = ppsc->dx[1] * fnqs / dt;
  creal fnqzs = ppsc->dx[2] * fnqs / dt;
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  fields_zero(pf, JXI);
  fields_zero(pf, JYI);
  fields_zero(pf, JZI);
  
  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    // x^n, p^n -> x^(n+.5), p^n

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
    int lg2 = fint(v);
    int lg3 = fint(w);
    creal h2 = v - lg2;
    creal h3 = w - lg3;

    creal g0y = 1.f - h2;
    creal g0z = 1.f - h3;
    creal g1y = h2;
    creal g1z = h3;

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0Y(+0) = g0y;
    S0Y(+1) = g1y;
    S0Z(+0) = g0z;
    S0Z(+1) = g1z;

    v = (part->yi - patch->xb[1]) * dyi - .5f;
    w = (part->zi - patch->xb[2]) * dzi - .5f;
    int lh2 = fint(v);
    int lh3 = fint(w);
    h2 = v - lh2;
    h3 = w - lh3;
    creal h0y = 1.f - h2;
    creal h0z = 1.f - h3;
    creal h1y = h2;
    creal h1z = h3;

    // FIELD INTERPOLATION

#define INTERPOLATE_FIELD(m, gy, gz)					\
    (gz##0z*(gy##0y*F3(m, 0,l##gy##2  ,l##gz##3  ) +			\
	     gy##1y*F3(m, 0,l##gy##2+1,l##gz##3  )) +			\
     gz##1z*(gy##0y*F3(m, 0,l##gy##2  ,l##gz##3+1) +			\
	     gy##1y*F3(m, 0,l##gy##2+1,l##gz##3+1)))			\
      
    creal exq = INTERPOLATE_FIELD(EX, g, g);
    creal eyq = INTERPOLATE_FIELD(EY, h, g);
    creal ezq = INTERPOLATE_FIELD(EZ, g, h);

    creal hxq = INTERPOLATE_FIELD(HX, h, h);
    creal hyq = INTERPOLATE_FIELD(HY, g, h);
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

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal yi = part->yi + vyi * yl;
    creal zi = part->zi + vzi * zl;

    v = (yi - patch->xb[1]) * dyi;
    w = (zi - patch->xb[2]) * dzi;
    int k2 = fint(v);
    int k3 = fint(w);
    h2 = v - k2;
    h3 = w - k3;

    for (int i = -1; i <= 2; i++) {
      S1Y(i) = 0.f;
      S1Z(i) = 0.f;
    }

    S1Y(k2-lg2+0) = 1.f - h2;
    S1Y(k2-lg2+1) = h2;
    S1Z(k3-lg3+0) = 1.f - h3;
    S1Z(k3-lg3+1) = h3;

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = 0; i <= 1; i++) {
      S1Y(i) -= S0Y(i);
      S1Z(i) -= S0Z(i);
    }

    int l2min, l3min, l2max, l3max;
    
    if (k2 == lg2) {
      l2min = 0; l2max = +1;
    } else if (k2 == lg2 - 1) {
      l2min = -1; l2max = +1;
    } else { // (k2 == lg2 + 1)
      l2min = 0; l2max = +2;
    }

    if (k3 == lg3) {
      l3min = 0; l3max = +1;
    } else if (k3 == lg3 - 1) {
      l3min = -1; l3max = +1;
    } else { // (k3 == lg3 + 1)
      l3min = 0; l3max = +2;
    }

    creal fnqx = vxi * part->qni * part->wni * fnqs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l2 = l2min; l2 <= l2max; l2++) {
	creal wx = S0Y(l2) * S0Z(l3)
	  + .5f * S1Y(l2) * S0Z(l3)
	  + .5f * S0Y(l2) * S1Z(l3)
	  + (1.f/3.f) * S1Y(l2) * S1Z(l3);
	creal jxh = fnqx*wx;
	F3(JXI, 0,lg2+l2,lg3+l3) += jxh;
      }
    }

    creal fnqy = part->qni * part->wni * fnqys;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      creal jyh = 0.f;
      for (int l2 = l2min; l2 < l2max; l2++) {
	creal wy = S1Y(l2) * (S0Z(l3) + .5f*S1Z(l3));
	jyh -= fnqy*wy;
	F3(JYI, 0,lg2+l2,lg3+l3) += jyh;
      }
    }

    creal fnqz = part->qni * part->wni * fnqzs;
    for (int l2 = l2min; l2 <= l2max; l2++) {
      creal jzh = 0.f;
      for (int l3 = l3min; l3 < l3max; l3++) {
	creal wz = S1Z(l3) * (S0Y(l2) + .5f*S1Y(l2));
	jzh -= fnqz*wz;
	F3(JZI, 0,lg2+l2,lg3+l3) += jzh;
      }
    }
  }
}

void
psc_push_particles_1st_push_yz(struct psc_push_particles *push,
			       mparticles_base_t *particles_base,
			       mfields_base_t *flds_base)
{
  mfields_t flds;
  mparticles_t particles;
  psc_mfields_get_from(&flds, EX, EX + 6, flds_base);
  psc_mparticles_get_from(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("1st_part_yz", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_push_part_1st_yz(p, &flds.f[p], &particles.p[p]);
  }
  prof_stop(pr);

  psc_mfields_put_to(&flds, JXI, JXI + 3, flds_base);
  psc_mparticles_put_to(&particles, particles_base);
}

