
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "inc_defs.h"

#define DIM DIM_YZ

#include "inc_params.c"
#include "inc_cache.c"
#include "inc_interpolate.c"
#include "inc_push.c"

static void
do_push_part_1st_yz(int p, fields_t *pf, particle_range_t prts)
{
#define S0Y(off) s0y[off+1]
#define S0Z(off) s0z[off+1]
#define S1Y(off) s1y[off+1]
#define S1Z(off) s1z[off+1]

  particle_real_t s0y[4] = {}, s0z[4] = {}, s1y[4], s1z[4];

  particle_real_t dt = ppsc->dt;
  particle_real_t dqs = .5f * ppsc->coeff.eta * dt;
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t fnqys = ppsc->patch[p].dx[1] * fnqs / dt;
  particle_real_t fnqzs = ppsc->patch[p].dx[2] * fnqs / dt;

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);
    particle_real_t vxi[3];

    // x^n, p^n -> x^(n+.5), p^n
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5f * dt);

    // field interpolation

    int lg[3], lh[3];
    particle_real_t og[3], oh[3];
    find_idx_off_1st_rel(&part->xi, lg, og, 0.f); // FIXME passing xi hack
    find_idx_off_1st_rel(&part->xi, lh, oh, -.5f);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0Y(+0) = 1.f - og[1];
    S0Y(+1) = og[1];
    S0Z(+0) = 1.f - og[2];
    S0Z(+1) = og[2];

    // FIELD INTERPOLATION
    particle_real_t exq, eyq, ezq, hxq, hyq, hzq;
    INTERPOLATE_1ST_STD(pf, exq, eyq, ezq, hxq, hyq, hzq);

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    particle_real_t dq = dqs * particle_qni_div_mni(part);
    push_pxi(part, exq, eyq, ezq, hxq, hyq, hzq, dq);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5 * dt);

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    particle_real_t xi[3] = { 0.f,
		    part->yi + vxi[1] * .5f * dt, 
		    part->zi + vxi[2] * .5f * dt };

    int lf[3];
    particle_real_t of[3];
    find_idx_off_1st_rel(xi, lf, of, 0.f);

    for (int i = -1; i <= 2; i++) {
      S1Y(i) = 0.f;
      S1Z(i) = 0.f;
    }

    S1Y(lf[1]-lg[1]+0) = 1.f - of[1];
    S1Y(lf[1]-lg[1]+1) = of[1];
    S1Z(lf[2]-lg[2]+0) = 1.f - of[2];
    S1Z(lf[2]-lg[2]+1) = of[2];

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = 0; i <= 1; i++) {
      S1Y(i) -= S0Y(i);
      S1Z(i) -= S0Z(i);
    }

    int l2min, l3min, l2max, l3max;
    
    if (lf[1] == lg[1]) {
      l2min = 0; l2max = +1;
    } else if (lf[1] == lg[1] - 1) {
      l2min = -1; l2max = +1;
    } else { // (lf[1] == lg[1] + 1)
      l2min = 0; l2max = +2;
    }

    if (lf[2] == lg[2]) {
      l3min = 0; l3max = +1;
    } else if (lf[2] == lg[2] - 1) {
      l3min = -1; l3max = +1;
    } else { // (lf[2] == lg[2] + 1)
      l3min = 0; l3max = +2;
    }

    particle_real_t fnqx = vxi[0] * part->qni * part->wni * fnqs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l2 = l2min; l2 <= l2max; l2++) {
	particle_real_t wx = S0Y(l2) * S0Z(l3)
	  + .5f * S1Y(l2) * S0Z(l3)
	  + .5f * S0Y(l2) * S1Z(l3)
	  + (1.f/3.f) * S1Y(l2) * S1Z(l3);
	particle_real_t jxh = fnqx*wx;
	F3(pf, JXI, 0,lg[1]+l2,lg[2]+l3) += jxh;
      }
    }

    particle_real_t fnqy = part->qni * part->wni * fnqys;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      particle_real_t jyh = 0.f;
      for (int l2 = l2min; l2 < l2max; l2++) {
	particle_real_t wy = S1Y(l2) * (S0Z(l3) + .5f*S1Z(l3));
	jyh -= fnqy*wy;
	F3(pf, JYI, 0,lg[1]+l2,lg[2]+l3) += jyh;
      }
    }

    particle_real_t fnqz = part->qni * part->wni * fnqzs;
    for (int l2 = l2min; l2 <= l2max; l2++) {
      particle_real_t jzh = 0.f;
      for (int l3 = l3min; l3 < l3max; l3++) {
	particle_real_t wz = S1Z(l3) * (S0Y(l2) + .5f*S1Y(l2));
	jzh -= fnqz*wz;
	F3(pf, JZI, 0,lg[1]+l2,lg[2]+l3) += jzh;
      }
    }
  }
}

void
psc_push_particles_1st_push_mprts_yz(struct psc_push_particles *push,
				     struct psc_mparticles *mprts,
				     struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("1st_part_yz", 1., 0, 0);
  }

  prof_start(pr);
  params_1vb_set(ppsc, NULL, NULL);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    particle_range_t prts = particle_range_mprts(mprts, p);

    psc_fields_zero_range(flds, JXI, JXI + 3);
    do_push_part_1st_yz(p, flds, prts);
  }
  prof_stop(pr);
}

