
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "c_common_push.c"

static void
do_push_part_1vb_yz(struct psc_fields *pf, struct psc_particles *pp)
{
  particle_real_t dt = ppsc->dt;
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t fnqys = ppsc->patch[pf->p].dx[1] * fnqs / dt;
  particle_real_t fnqzs = ppsc->patch[pf->p].dx[2] * fnqs / dt;
  particle_real_t dxi[3] = { 1.f / ppsc->patch[pf->p].dx[0], 1.f / ppsc->patch[pf->p].dx[1], 1.f / ppsc->patch[pf->p].dx[2] };
  particle_real_t dq_kind[ppsc->nr_kinds];
  particle_real_t fnqx_kind[ppsc->nr_kinds];
  particle_real_t fnqy_kind[ppsc->nr_kinds];
  particle_real_t fnqz_kind[ppsc->nr_kinds];
  for (int k = 0; k < ppsc->nr_kinds; k++) {
    dq_kind[k] = .5f * ppsc->coeff.eta * dt * ppsc->kinds[k].q / ppsc->kinds[k].m;
    fnqx_kind[k] = fnqs * ppsc->kinds[k].q;
    fnqy_kind[k] = fnqys * ppsc->kinds[k].q;
    fnqz_kind[k] = fnqzs * ppsc->kinds[k].q;
  }

  struct psc_particles_single *sngl = psc_particles_single(pp);
  memset(sngl->b_cnt, 0, (sngl->nr_blocks + 1) * sizeof(*sngl->b_cnt));

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    // field interpolation

    int lg[3], lh[3];
    particle_real_t og[3], oh[3], xm[3];
    find_idx_off_pos_1st_rel(&part->xi, lg, og, xm, 0.f, dxi); // FIXME passing xi hack
    find_idx_off_1st_rel(&part->xi, lh, oh, -.5f, dxi);

    // FIELD INTERPOLATION
    particle_real_t exq, eyq, ezq, hxq, hyq, hzq;
    INTERPOLATE_1ST_STD(pf, exq, eyq, ezq, hxq, hyq, hzq);

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    particle_real_t dq = dq_kind[part->kind];
    push_pxi(part, exq, eyq, ezq, hxq, hyq, hzq, dq);

    particle_real_t vxi[3];

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5f * dt);

    // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt

    int lf[3];
    particle_real_t of[3];
    find_idx_off_1st_rel(&part->xi, lf, of, 0.f, dxi);

    particle_real_t fnqx = vxi[0] * particle_wni(part) * fnqx_kind[part->kind];
    F3_CURR(pf, JXI, 0,lf[1]  ,lf[2]  ) += (1.f - of[1]) * (1.f - of[2]) * fnqx;
    F3_CURR(pf, JXI, 0,lf[1]+1,lf[2]  ) += (      of[1]) * (1.f - of[2]) * fnqx;
    F3_CURR(pf, JXI, 0,lf[1]  ,lf[2]+1) += (1.f - of[1]) * (      of[2]) * fnqx;
    F3_CURR(pf, JXI, 0,lf[1]+1,lf[2]+1) += (      of[1]) * (      of[2]) * fnqx;

    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5f * dt);

    particle_real_t xp[3];
    find_idx_off_pos_1st_rel(&part->xi, lf, of, xp, 0.f, dxi);

    // FIXME, only if blocksize == 1!
    int *b_pos = lf;
    int *b_mx = sngl->b_mx;
    if (b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      sngl->b_idx[n] = b_pos[2] * b_mx[1] + b_pos[1];
    } else { // out of bounds
      sngl->b_idx[n] = sngl->nr_blocks;
      assert(sngl->b_cnt[sngl->nr_blocks] < sngl->n_alloced);
      // append to back
      *particles_get_one(pp, pp->n_part + sngl->b_cnt[sngl->nr_blocks]) = *part;
    }
    sngl->b_cnt[sngl->b_idx[n]]++;

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

    particle_real_t fnq[2] = { particle_wni(part) * fnqy_kind[part->kind],
			       particle_wni(part) * fnqz_kind[part->kind] };

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

static void
psc_push_particles_push_a_yz(struct psc_push_particles *push,
			     struct psc_particles *prts_base,
			     struct psc_fields *flds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register(PARTICLE_TYPE "_1vb_push_yz", 1., 0, 0);
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

