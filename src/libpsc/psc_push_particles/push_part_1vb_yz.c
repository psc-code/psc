
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "c_common_push.c"

static inline void
calc_dx1(creal dx1[2], creal x[2], creal dx[2], int off[2])
{
  if (off[1] == 0) {
    dx1[0] = .5 * off[0] - x[0];
    dx1[1] = dx[1] / dx[0] * dx1[0];
  } else {
    dx1[1] = .5 * off[1] - x[1];
    dx1[0] = dx[0] / dx[1] * dx1[1];
  }
}

static inline void
curr_2d_vb_cell(fields_t *pf, int i[2], creal x[2], creal dx[2], creal fnq[2],
		creal dxt[2], int off[2])
{
  F3(pf, JYI, 0,i[0],i[1]  ) += fnq[0] * dx[0] * (.5 - x[1] - .5 * dx[1]);
  F3(pf, JYI, 0,i[0],i[1]+1) += fnq[0] * dx[0] * (.5 + x[1] + .5 * dx[1]);
  F3(pf, JZI, 0,i[0],i[1]  ) += fnq[1] * dx[1] * (.5 - x[0] - .5 * dx[0]);
  F3(pf, JZI, 0,i[0]+1,i[1]) += fnq[1] * dx[1] * (.5 + x[0] + .5 * dx[0]);
  if (dxt) {
    dxt[0] -= dx[0];
    dxt[1] -= dx[1];
    x[0] += dx[0] - off[0];
    x[1] += dx[1] - off[1];
    i[0] += off[0];
    i[1] += off[1];
  }
}

static void
do_push_part_1vb_yz(int p, fields_t *pf, particles_t *pp)
{
  creal dt = ppsc->dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqys = ppsc->dx[1] * fnqs / dt;
  creal fnqzs = ppsc->dx[2] * fnqs / dt;
  creal dxi[3] = { 1.f / ppsc->dx[0], 1.f / ppsc->dx[1], 1.f / ppsc->dx[2] };

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
    creal vxi[3];

    // x^n, p^n -> x^(n+.5), p^n
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5f * dt);

    // field interpolation

    int lg[3], lh[3];
    creal og[3], oh[3], xm[3];
    find_idx_off_pos_1st(&part->xi, lg, og, xm, 0.f, patch->xb, dxi); // FIXME passing xi hack
    find_idx_off_1st(&part->xi, lh, oh, -.5f, patch->xb, dxi);

    // FIELD INTERPOLATION

    INTERPOLATE_SETUP_1ST;

    creal exq = INTERPOLATE_FIELD_1ST(EX, g, g);
    creal eyq = INTERPOLATE_FIELD_1ST(EY, h, g);
    creal ezq = INTERPOLATE_FIELD_1ST(EZ, g, h);

    creal hxq = INTERPOLATE_FIELD_1ST(HX, h, h);
    creal hyq = INTERPOLATE_FIELD_1ST(HY, g, h);
    creal hzq = INTERPOLATE_FIELD_1ST(HZ, h, g);

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    push_pxi(part, exq, eyq, ezq, hxq, hyq, hzq, dqs);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5 * dt);

    // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt

    int lf[3];
    creal of[3];
    find_idx_off_1st(&part->xi, lf, of, 0.f, patch->xb, dxi);

    creal fnqx = vxi[0] * part->qni * part->wni * fnqs;
    F3(pf, JXI, 0,lf[1]  ,lf[2]  ) += (1.f - of[1]) * (1.f - of[2]) * fnqx;
    F3(pf, JXI, 0,lf[1]+1,lf[2]  ) += (      of[1]) * (1.f - of[2]) * fnqx;
    F3(pf, JXI, 0,lf[1]  ,lf[2]+1) += (1.f - of[1]) * (      of[2]) * fnqx;
    F3(pf, JXI, 0,lf[1]+1,lf[2]+1) += (      of[1]) * (      of[2]) * fnqx;

    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal xi[3] = { 0.f,
		    part->yi + vxi[1] * .5f * dt,
		    part->zi + vxi[2] * .5f * dt };

    creal xp[3];
    find_idx_off_pos_1st(xi, lf, of, xp, 0.f, patch->xb, dxi);

    // OUT OF PLANE CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt

    int i[2] = { lg[1], lg[2] };
    int idiff[2] = { lf[1] - lg[1], lf[2] - lg[2] };
    creal dx[2] = { xp[1] - xm[1], xp[2] - xm[2] };
    creal x[2] = { xm[1] - (i[0] + .5), xm[2] - (i[1] + .5) };

    creal dx1[2];
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
      dx1[0] = .5 * idiff[0] - x[0];
      dx1[1] = dx[1] / dx[0] * dx1[0];
      if (creal_abs(x[1] + dx1[1]) > .5f) {
	first_dir = 1;
      } else {
	first_dir = 0;
      }
      second_dir = 1 - first_dir;
    }

    creal fnq[2] = { part->qni * part->wni * fnqys,
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
psc_push_particles_1vb_push_yz(struct psc_push_particles *push,
			       mparticles_base_t *particles_base,
			       mfields_base_t *flds_base)
{
  mparticles_t particles;
  psc_mparticles_get_from(&particles, particles_base);
  mfields_t *flds = psc_mfields_get_from(EX, EX + 6, flds_base);

  static int pr;
  if (!pr) {
    pr = prof_register("1vb_push_yz", 1., 0, 0);
  }
  prof_start(pr);
  psc_mfields_zero(flds, JXI);
  psc_mfields_zero(flds, JYI);
  psc_mfields_zero(flds, JZI);

  psc_foreach_patch(ppsc, p) {
    do_push_part_1vb_yz(p, &flds->f[p], &particles.p[p]);
  }
  prof_stop(pr);

  psc_mfields_put_to(flds, JXI, JXI + 3, flds_base);
  psc_mparticles_put_to(&particles, particles_base);
}

