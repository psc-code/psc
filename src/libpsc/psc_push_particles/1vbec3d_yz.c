
#include "psc_debug.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PUSHER_TYPE "1vbec3d"

#include "c_common_push.c"

static void
do_push_part_1vb_yz(struct psc_fields *pf, struct psc_particles *pp)
{
  particle_real_t dt = ppsc->dt;
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t fnqxs = ppsc->patch[pf->p].dx[0] * fnqs / dt;
  particle_real_t fnqys = ppsc->patch[pf->p].dx[1] * fnqs / dt;
  particle_real_t fnqzs = ppsc->patch[pf->p].dx[2] * fnqs / dt;
  particle_real_t dxi[3] = { 1.f / ppsc->patch[pf->p].dx[0], 1.f / ppsc->patch[pf->p].dx[1], 1.f / ppsc->patch[pf->p].dx[2] };
  particle_real_t dq_kind[ppsc->nr_kinds];
  particle_real_t fnqx_kind[ppsc->nr_kinds];
  particle_real_t fnqy_kind[ppsc->nr_kinds];
  particle_real_t fnqz_kind[ppsc->nr_kinds];
  for (int k = 0; k < ppsc->nr_kinds; k++) {
    dq_kind[k] = .5f * ppsc->coeff.eta * dt * ppsc->kinds[k].q / ppsc->kinds[k].m;
    fnqx_kind[k] = fnqxs * ppsc->kinds[k].q;
    fnqy_kind[k] = fnqys * ppsc->kinds[k].q;
    fnqz_kind[k] = fnqzs * ppsc->kinds[k].q;
  }

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    // field interpolation

    int lg[3], lh[3];
    particle_real_t og[3], oh[3], xm[3];
    find_idx_off_pos_1st_rel(&part->xi, lg, og, xm, 0.f, dxi); // FIXME passing xi hack
    find_idx_off_1st_rel(&part->xi, lh, oh, -.5f, dxi);

    // FIELD INTERPOLATION
    particle_real_t exq, eyq, ezq, hxq, hyq, hzq;
    INTERPOLATE_1ST_EC(exq, eyq, ezq, hxq, hyq, hzq);

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
    particle_real_t dq = dq_kind[part->kind];
    push_pxi(part, exq, eyq, ezq, hxq, hyq, hzq, dq);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
    particle_real_t vxi[3];
    calc_vxi(vxi, part);
    push_xi(part, vxi, dt);

    int lf[3];
    particle_real_t of[3], xp[3];
    find_idx_off_pos_1st_rel(&part->xi, lf, of, xp, 0.f, dxi);

    // IN PLANE CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt

    int i[3] = { 0, lg[1], lg[2] };
    int idiff[3] = { 0, lf[1] - lg[1], lf[2] - lg[2] };
    particle_real_t dx[3] = { 0., xp[1] - xm[1], xp[2] - xm[2] };
    particle_real_t x[3] = { 0., xm[1] - (i[1] + .5f), xm[2] - (i[2] + .5f) };

    particle_real_t dx1[3];
    int off[3];
    int first_dir, second_dir = -1;
    // FIXME, make sure we never div-by-zero?
    if (idiff[1] == 0 && idiff[2] == 0) {
      first_dir = -1;
    } else if (idiff[1] == 0) {
      first_dir = 2;
    } else if (idiff[2] == 0) {
      first_dir = 1;
    } else {
      dx1[1] = .5f * idiff[1] - x[1];
     if (dx[1] == 0.f) {
      dx1[2] = 0.f;
     } else {
      dx1[2] = dx[2] / dx[1] * dx1[1];
     }
      if (particle_real_abs(x[2] + dx1[2]) > .5f) {
	first_dir = 2;
      } else {
	first_dir = 1;
      }
      second_dir = 3 - first_dir;
    }

    particle_real_t fnq[3] = { particle_wni(part) * fnqx_kind[part->kind],
			       particle_wni(part) * fnqy_kind[part->kind],
			       particle_wni(part) * fnqz_kind[part->kind] };
    dx[0] = vxi[0] * dt * dxi[0];

    if (first_dir >= 0) {
      off[3 - first_dir] = 0;
      off[first_dir] = idiff[first_dir];
      calc_3d_dx1(dx1, x, dx, off);
      curr_3d_vb_cell(pf, i, x, dx1, fnq, dx, off);
    }

    if (second_dir >= 0) {
      off[first_dir] = 0;
      off[second_dir] = idiff[second_dir];
      calc_3d_dx1(dx1, x, dx, off);
      curr_3d_vb_cell(pf, i, x, dx1, fnq, dx, off);
    }
    
    curr_3d_vb_cell(pf, i, x, dx, fnq, NULL, NULL);
  }
}

void
psc_push_particles_push_yz(struct psc_push_particles *push,
			   struct psc_particles *prts_base,
			   struct psc_fields *flds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register(PARTICLE_TYPE "_" PUSHER_TYPE "_push_yz", 1., 0, 0);
  }

  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 6);
  
  prof_start(pr);
  psc_fields_zero_range(flds, JXI, JXI + 3);
  struct psc_fields *flds_cache = cache_fields_from_em(flds);
  do_push_part_1vb_yz(flds_cache, prts);
  cache_fields_to_j(flds_cache, flds);
  psc_fields_destroy(flds_cache);
  prof_stop(pr);

  psc_particles_put_as(prts, prts_base, 0);
  psc_fields_put_as(flds, flds_base, JXI, JXI + 3);
}

