
#include "psc_push_particles_1st.h"
#include "psc_debug.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PARTICLE_TYPE "c"

#include "inc_defs.h"

#define DIM DIM_YZ
#define CALC_J CALC_J_1VB_VAR1
#define F3_CURR F3_C
#define F3_CACHE F3_C
#define F3_CACHE_TYPE "c"

#include "inc_params.c"
#include "inc_interpolate.c"
#include "inc_push.c"
#include "inc_curr.c"
#include "c_common_push.c"

static void
do_push_part_1vb_yz(struct psc_fields *pf, struct psc_particles *pp)
{
  int p = pf->p;
  particle_real_t dt = ppsc->dt;
  particle_real_t dqs = .5f * ppsc->coeff.eta * dt;
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t fnqys = ppsc->patch[p].dx[1] * fnqs / dt;
  particle_real_t fnqzs = ppsc->patch[p].dx[2] * fnqs / dt;
  particle_real_t dxi[3] = { 1.f / ppsc->patch[p].dx[0], 1.f / ppsc->patch[p].dx[1], 1.f / ppsc->patch[p].dx[2] };

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
    particle_real_t vxi[3];

    // x^n, p^n -> x^(n+.5), p^n
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5f * dt);

    // field interpolation

    int lg[3], lh[3];
    particle_real_t og[3], oh[3], xm[3];
    find_idx_off_pos_1st_rel(&part->xi, lg, og, xm, 0.f, dxi); // FIXME passing xi hack
    find_idx_off_1st_rel(&part->xi, lh, oh, -.5f, dxi);

    // FIELD INTERPOLATION

    particle_real_t exq, eyq, ezq, hxq, hyq, hzq;
    INTERPOLATE_1ST_STD(pf, exq, eyq, ezq, hxq, hyq, hzq);

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    particle_real_t dq = dqs * particle_qni_div_mni(part);
    push_pxi(part, exq, eyq, ezq, hxq, hyq, hzq, dq);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5f * dt);

    // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt

    int lf[3];
    particle_real_t of[3];
    find_idx_off_1st_rel(&part->xi, lf, of, 0.f, dxi);

    particle_real_t fnqx = vxi[0] * part->qni * part->wni * fnqs;
    F3_CURR(pf, JXI, 0,lf[1]  ,lf[2]  ) += (1.f - of[1]) * (1.f - of[2]) * fnqx;
    F3_CURR(pf, JXI, 0,lf[1]+1,lf[2]  ) += (      of[1]) * (1.f - of[2]) * fnqx;
    F3_CURR(pf, JXI, 0,lf[1]  ,lf[2]+1) += (1.f - of[1]) * (      of[2]) * fnqx;
    F3_CURR(pf, JXI, 0,lf[1]+1,lf[2]+1) += (      of[1]) * (      of[2]) * fnqx;

    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    particle_real_t xi[3] = { 0.f,
		    part->yi + vxi[1] * .5f * dt,
		    part->zi + vxi[2] * .5f * dt };

    particle_real_t xp[3];
    find_idx_off_pos_1st_rel(xi, lf, of, xp, 0.f, dxi);

    // OUT OF PLANE CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt

    int i[2] = { lg[1], lg[2] };
    int idiff[2] = { lf[1] - lg[1], lf[2] - lg[2] };
    particle_real_t dx[2] = { xp[1] - xm[1], xp[2] - xm[2] };
    particle_real_t x[2] = { xm[1] - (i[0] + .5f), xm[2] - (i[1] + .5f) };

    particle_real_t dx1[2];
    int off[2];
    int first_dir, second_dir = -1;
    // FIXME, make sure we never div-by-zero?
    if (idiff[0] == 0 && idiff[1] == 0) {
      first_dir = -1;
    } else if (idiff[0] == 0) {
      first_dir = 1;
    } else if (idiff[1] == 0) {
      first_dir = 0;
    } else {
      dx1[0] = .5f * idiff[0] - x[0];
      dx1[1] = dx[1] / dx[0] * dx1[0];
      if (particle_real_abs(x[1] + dx1[1]) > .5f) {
	first_dir = 1;
      } else {
	first_dir = 0;
      }
      second_dir = 1 - first_dir;
    }

    particle_real_t fnq[2] = { part->qni * part->wni * fnqys,
			       part->qni * part->wni * fnqzs };

    if (first_dir >= 0) {
      off[1-first_dir] = 0;
      off[first_dir] = idiff[first_dir];
      calc_dx1(dx1, x, dx, off);
      curr_2d_vb_cell(pf, i, x, dx1, fnq, dx, off);
    }

    if (second_dir >= 0) {
      off[first_dir] = 0;
      off[second_dir] = idiff[second_dir];
      calc_dx1(dx1, x, dx, off);
      curr_2d_vb_cell(pf, i, x, dx1, fnq, dx, off);
    }
    
    curr_2d_vb_cell(pf, i, x, dx, fnq, NULL, NULL);
  }
}

void
psc_push_particles_1vb_c_push_a_yz(struct psc_push_particles *push,
				   struct psc_particles *prts_base,
				   struct psc_fields *flds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("1vb_" PARTICLE_TYPE "_push_yz", 1., 0, 0);
  }

  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 6);
  
  prof_start(pr);
  params_1vb_set(ppsc, flds->p);
  psc_fields_zero_range(flds, JXI, JXI + 3);
  struct psc_fields *flds_cache = cache_fields_from_em(flds);
  do_push_part_1vb_yz(flds_cache, prts);
  cache_fields_to_j(flds_cache, flds);
  psc_fields_destroy(flds_cache);
  prof_stop(pr);
  
  psc_particles_put_as(prts, prts_base, 0);
  psc_fields_put_as(flds, flds_base, JXI, JXI + 3);
}

