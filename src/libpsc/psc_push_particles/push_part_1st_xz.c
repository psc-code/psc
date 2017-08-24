
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIM DIM_XZ
#define ORDER ORDER_1ST
#include "push_part_common.c"

static void
do_push_part_1st_xz(int p, struct psc_fields *pf, particle_range_t prts)
{
#include "push_part_common_vars.c"

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);

    // x^n, p^n -> x^(n+.5), p^n

    particle_real_t vv[3];
    calc_v(vv, &part->pxi);
    push_x(&part->xi, vv);

    particle_real_t u = part->xi * dxi;
    particle_real_t w = part->zi * dzi;
    int lg1 = particle_real_fint(u);
    int lg3 = particle_real_fint(w);
    particle_real_t h1 = u - lg1;
    particle_real_t h3 = w - lg3;

    particle_real_t g0x = 1.f - h1;
    particle_real_t g0z = 1.f - h3;
    particle_real_t g1x = h1;
    particle_real_t g1z = h3;

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S(s0x, 0) = g0x;
    S(s0x, 1) = g1x;
    S(s0z, 0) = g0z;
    S(s0z, 1) = g1z;

    u = part->xi * dxi - .5f;
    w = part->zi * dzi - .5f;
    int lh1 = particle_real_fint(u);
    int lh3 = particle_real_fint(w);
    h1 = u - lh1;
    h3 = w - lh3;
    particle_real_t h0x = 1.f - h1;
    particle_real_t h0z = 1.f - h3;
    particle_real_t h1x = h1;
    particle_real_t h1z = h3;

    // FIELD INTERPOLATION

#define INTERPOLATE_FIELD(m, gx, gz)					\
    (gz##0z*(gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3  ) +			\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3  )) +			\
     gz##1z*(gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3+1) +			\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3+1)))			\
      
    particle_real_t exq = INTERPOLATE_FIELD(EX, h, g);
    particle_real_t eyq = INTERPOLATE_FIELD(EY, g, g);
    particle_real_t ezq = INTERPOLATE_FIELD(EZ, g, h);

    particle_real_t hxq = INTERPOLATE_FIELD(HX, g, h);
    particle_real_t hyq = INTERPOLATE_FIELD(HY, h, h);
    particle_real_t hzq = INTERPOLATE_FIELD(HZ, h, g);

     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    particle_real_t dq = dqs * part->qni / part->mni;
    particle_real_t pxm = part->pxi + dq*exq;
    particle_real_t pym = part->pyi + dq*eyq;
    particle_real_t pzm = part->pzi + dq*ezq;

    particle_real_t root = dq / particle_real_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
    particle_real_t taux = hxq*root;
    particle_real_t tauy = hyq*root;
    particle_real_t tauz = hzq*root;

    particle_real_t tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    particle_real_t pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.f*taux*tauy+2.f*tauz)*pym + 
		(2.f*taux*tauz-2.f*tauy)*pzm)*tau;
    particle_real_t pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
		(1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.f*tauy*tauz+2.f*taux)*pzm)*tau;
    particle_real_t pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
		(2.f*tauy*tauz-2.f*taux)*pym +
		(1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
    
    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    calc_v(vv, &part->pxi);
    push_x(&part->xi, vv);

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    particle_real_t xi = part->xi + vv[0] * xl;
    particle_real_t zi = part->zi + vv[2] * zl;

    u = xi * dxi;
    w = zi * dzi;
    int k1 = particle_real_fint(u);
    int k3 = particle_real_fint(w);
    h1 = u - k1;
    h3 = w - k3;

    for (int i = -1; i <= 2; i++) {
      S(s1x, i) = 0.f;
      S(s1z, i) = 0.f;
    }

    S(s1x, k1-lg1+0) = 1.f - h1;
    S(s1x, k1-lg1+1) = h1;
    S(s1z, k3-lg3+0) = 1.f - h3;
    S(s1z, k3-lg3+1) = h3;

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = 0; i <= 1; i++) {
      S(s1x, i) -= S(s0x, i);
      S(s1z, i) -= S(s0z, i);
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

    particle_real_t fnqx = part->qni * part->wni * fnqxs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      particle_real_t jxh = 0.f;
      for (int l1 = l1min; l1 < l1max; l1++) {
	particle_real_t wx = S(s1x, l1) * (S(s0z, l3) + .5f*S(s1z, l3));
	jxh -= fnqx*wx;
	F3(pf, JXI, lg1+l1,0,lg3+l3) += jxh;
      }
    }

    particle_real_t fnqy = vv[1] * part->qni * part->wni * fnqs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l1 = l1min; l1 <= l1max; l1++) {
	particle_real_t wy = S(s0x, l1) * S(s0z, l3)
	  + .5f * S(s1x, l1) * S(s0z, l3)
	  + .5f * S(s0x, l1) * S(s1z, l3)
	  + (1.f/3.f) * S(s1x, l1) * S(s1z, l3);
	particle_real_t jyh = fnqy*wy;
	F3(pf, JYI, lg1+l1,0,lg3+l3) += jyh;
      }
    }

    particle_real_t fnqz = part->qni * part->wni * fnqzs;
    for (int l1 = l1min; l1 <= l1max; l1++) {
      particle_real_t jzh = 0.f;
      for (int l3 = l3min; l3 < l3max; l3++) {
	particle_real_t wz = S(s1z, l3) * (S(s0x, l1) + .5f*S(s1x, l1));
	jzh -= fnqz*wz;
	F3(pf, JZI, lg1+l1,0,lg3+l3) += jzh;
      }
    }
  }
}

void
psc_push_particles_1st_push_mprts_xz(struct psc_push_particles *push,
				     struct psc_mparticles *mprts,
				     struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register(PARTICLE_TYPE "_1st_push_xz", 1., 0, 0);
  }

  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    particle_range_t prts = particle_range_mprts(mprts, p);
    psc_fields_zero_range(flds, JXI, JXI + 3);
    do_push_part_1st_xz(p, flds, prts);
  }
  prof_stop(pr);
}

