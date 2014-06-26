
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "c_common_push.c"

static inline void
calc_3d_dx1(particle_real_t dx1[3], particle_real_t x[3], particle_real_t dx[3], int off[3])
{
  if (off[2] == 0) {
    dx1[1] = .5f * off[1] - x[1];
   if (dx[1] == 0.f) {
     dx1[0] = 0.f;
     dx1[2] = 0.f;
   } else {
     dx1[0] = dx[0] / dx[1] * dx1[1];
     dx1[2] = dx[2] / dx[1] * dx1[1];
   }
  } else {
    dx1[2] = .5f * off[2] - x[2];
   if (dx[2] == 0.f) {
     dx1[0] = 0.f;
     dx1[1] = 0.f;
   } else {
     dx1[0] = dx[0] / dx[2] * dx1[2];
     dx1[1] = dx[1] / dx[2] * dx1[2];
   }
  }
}

static inline void
curr_3d_vb_cell(struct psc_fields *pf, int i[3], particle_real_t x[3], particle_real_t dx[3],
		particle_real_t fnq[3], particle_real_t dxt[3], int off[3])
{
  particle_real_t h = (1.f/12.f) * dx[0] * dx[1] * dx[2];
  particle_real_t xa[3] = { 0.,
			    x[1] + .5f * dx[1],
			    x[2] + .5f * dx[2], };
  F3_CURR(pf, JXI, 0,i[1]  ,i[2]  ) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f - xa[2]) + h);
  F3_CURR(pf, JXI, 0,i[1]+1,i[2]  ) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f - xa[2]) - h);
  F3_CURR(pf, JXI, 0,i[1]  ,i[2]+1) += fnq[0] * (dx[0] * (.5f - xa[1]) * (.5f + xa[2]) + h);
  F3_CURR(pf, JXI, 0,i[1]+1,i[2]+1) += fnq[0] * (dx[0] * (.5f + xa[1]) * (.5f + xa[2]) - h);

  F3_CURR(pf, JYI, 0,i[1]  ,i[2]  ) += fnq[1] * dx[1] * (.5f - xa[2]);
  F3_CURR(pf, JYI, 0,i[1]  ,i[2]+1) += fnq[1] * dx[1] * (.5f + xa[2]);
  F3_CURR(pf, JZI, 0,i[1]  ,i[2]  ) += fnq[2] * dx[2] * (.5f - xa[1]);
  F3_CURR(pf, JZI, 0,i[1]+1,i[2]  ) += fnq[2] * dx[2] * (.5f + xa[1]);
  if (dxt) {
    dxt[0] -= dx[0];
    dxt[1] -= dx[1];
    dxt[2] -= dx[2];
    x[1] += dx[1] - off[1];
    x[2] += dx[2] - off[2];
    i[1] += off[1];
    i[2] += off[2];
  }
}

static void
do_push_part_1vbec_yz(struct psc_fields *pf, struct psc_particles *pp)
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

    INTERPOLATE_SETUP_1ST_EC;

    particle_real_t exq = 
      (g0z*(g0y*F3_CACHE(pf, EX, 0,lg[1]  ,lg[2]  ) +
	    g1y*F3_CACHE(pf, EX, 0,lg[1]+1,lg[2]  )) +
       g1z*(g0y*F3_CACHE(pf, EX, 0,lg[1]  ,lg[2]+1) +
	    g1y*F3_CACHE(pf, EX, 0,lg[1]+1,lg[2]+1)));
    particle_real_t eyq =
      (g0z*F3_CACHE(pf, EY, 0,lg[1]  ,lg[2]  ) +
       g1z*F3_CACHE(pf, EY, 0,lg[1]  ,lg[2]+1));
    particle_real_t ezq =
      (g0y*F3_CACHE(pf, EZ, 0,lg[1]  ,lg[2]  ) +
       g1y*F3_CACHE(pf, EZ, 0,lg[1]+1,lg[2]  ));

    particle_real_t hxq =
      F3_CACHE(pf, HX, 0,lg[1]  ,lg[2]  );
    particle_real_t hyq =
      (g0y*F3_CACHE(pf, HY, 0,lg[1]  ,lg[2]  ) +
       g1y*F3_CACHE(pf, HY, 0,lg[1]+1,lg[2]  ));
    particle_real_t hzq =
      (g0z*F3_CACHE(pf, HZ, 0,lg[1]  ,lg[2]  ) +
       g1z*F3_CACHE(pf, HZ, 0,lg[1]  ,lg[2]+1));

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
    particle_real_t dq = dq_kind[part->kind];
    push_pxi(part, exq, eyq, ezq, hxq, hyq, hzq, dq);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
    particle_real_t vxi[3];
    calc_vxi(vxi, part);
    push_xi(part, vxi, dt);

    particle_real_t xp[3];
    int lf[3];
    particle_real_t of[3];
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
psc_push_particles_1vbec_push_a_yz(struct psc_push_particles *push,
				   struct psc_particles *prts_base,
				   struct psc_fields *flds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register(PARTICLE_TYPE "_1vbec_push_yz", 1., 0, 0);
  }

  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 6);
  
  prof_start(pr);
  psc_fields_zero_range(flds, JXI, JXI + 3);
  struct psc_fields *flds_cache = cache_fields_from_em(flds);
  do_push_part_1vbec_yz(flds_cache, prts);
  cache_fields_to_j(flds_cache, flds);
  psc_fields_destroy(flds_cache);
  prof_stop(pr);

  psc_particles_put_as(prts, prts_base, 0);
  psc_fields_put_as(flds, flds_base, JXI, JXI + 3);
}
