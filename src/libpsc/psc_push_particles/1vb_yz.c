
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

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    // field interpolation

    int lg[3], lh[3];
    particle_real_t og[3], oh[3], xm[3];
    find_idx_off_pos_1st_rel(&part->xi, lg, og, xm, 0.f, dxi); // FIXME passing xi hack
    find_idx_off_1st_rel(&part->xi, lh, oh, -.5f, dxi);

    // FIELD INTERPOLATION
    particle_real_t exq, eyq, ezq, hxq, hyq, hzq;
    INTERPOLATE_1ST(exq, eyq, ezq, hxq, hyq, hzq);

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
    particle_real_t dq = dq_kind[part->kind];
    push_pxi(part, exq, eyq, ezq, hxq, hyq, hzq, dq);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
    particle_real_t vxi[3];
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5f * dt);

    // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt
    CALC_JX_2D(pf, part, vxi);

    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)
    push_xi(part, vxi, .5f * dt);

    int lf[3];
    particle_real_t of[3], xp[3];
    find_idx_off_pos_1st_rel(&part->xi, lf, of, xp, 0.f, dxi);

    // IN PLANE CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
    CALC_JYZ_2D(pf, xm, xp);
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

