
#include "psc_push_particles_single.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef fields_single_t fields_curr_t;
#define F3_CURR F3_S

#include "c_common_push.c"

static void
do_push_part_1vb_yz(int p, fields_single_t *pf, particles_t *pp)
{
  creal dt = ppsc->dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqys = ppsc->dx[1] * fnqs / dt;
  creal fnqzs = ppsc->dx[2] * fnqs / dt;
  creal dxi[3] = { 1.f / ppsc->dx[0], 1.f / ppsc->dx[1], 1.f / ppsc->dx[2] };

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
    creal vxi[3];

    // x^n, p^n -> x^(n+.5), p^n
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5f * dt);

    // field interpolation

    int lg[3], lh[3];
    creal og[3], oh[3], xm[3];
    find_idx_off_pos_1st_rel(&part->xi, lg, og, xm, 0.f, dxi); // FIXME passing xi hack
    find_idx_off_1st_rel(&part->xi, lh, oh, -.5f, dxi);

    // FIELD INTERPOLATION

    INTERPOLATE_SETUP_1ST;

#define INTERPOLATE_FIELD_1ST_C(m, gy, gz)				\
    (gz##0z*(gy##0y*F3_C(pf, m, 0,l##gy[1]  ,l##gz[2]  ) +		\
	     gy##1y*F3_C(pf, m, 0,l##gy[1]+1,l##gz[2]  )) +		\
     gz##1z*(gy##0y*F3_C(pf, m, 0,l##gy[1]  ,l##gz[2]+1) +		\
	     gy##1y*F3_C(pf, m, 0,l##gy[1]+1,l##gz[2]+1)))

    creal exq = INTERPOLATE_FIELD_1ST_C(EX, g, g);
    creal eyq = INTERPOLATE_FIELD_1ST_C(EY, h, g);
    creal ezq = INTERPOLATE_FIELD_1ST_C(EZ, g, h);

    creal hxq = INTERPOLATE_FIELD_1ST_C(HX, h, h);
    creal hyq = INTERPOLATE_FIELD_1ST_C(HY, g, h);
    creal hzq = INTERPOLATE_FIELD_1ST_C(HZ, h, g);

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    push_pxi(part, exq, eyq, ezq, hxq, hyq, hzq, dqs);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    calc_vxi(vxi, part);
    push_xi(part, vxi, .5 * dt);

    // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt

    int lf[3];
    creal of[3];
    find_idx_off_1st_rel(&part->xi, lf, of, 0.f, dxi);

    creal fnqx = vxi[0] * part->qni * part->wni * fnqs;
    F3_S(pf, JXI, 0,lf[1]  ,lf[2]  ) += (1.f - of[1]) * (1.f - of[2]) * fnqx;
    F3_S(pf, JXI, 0,lf[1]+1,lf[2]  ) += (      of[1]) * (1.f - of[2]) * fnqx;
    F3_S(pf, JXI, 0,lf[1]  ,lf[2]+1) += (1.f - of[1]) * (      of[2]) * fnqx;
    F3_S(pf, JXI, 0,lf[1]+1,lf[2]+1) += (      of[1]) * (      of[2]) * fnqx;

    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal xi[3] = { 0.f,
		    part->yi + vxi[1] * .5f * dt,
		    part->zi + vxi[2] * .5f * dt };

    creal xp[3];
    find_idx_off_pos_1st_rel(xi, lf, of, xp, 0.f, dxi);

    // OUT OF PLANE CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt

    int i[2] = { lg[1], lg[2] };
    int idiff[2] = { lf[1] - lg[1], lf[2] - lg[2] };
    creal dx[2] = { xp[1] - xm[1], xp[2] - xm[2] };
    creal x[2] = { xm[1] - (i[0] + .5f), xm[2] - (i[1] + .5f) };

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
      dx1[0] = .5f * idiff[0] - x[0];
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
psc_push_particles_single_1vb_push_yz(struct psc_push_particles *push,
				      mparticles_base_t *particles_base,
				      mfields_base_t *flds_base)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);
  mfields_t *flds = psc_mfields_get_cf(flds_base, EX, EX + 6);

  static int pr;
  if (!pr) {
    pr = prof_register("single_1vb_push_yz", 1., 0, 0);
  }
  prof_start(pr);
  psc_mfields_zero(flds, JXI);
  psc_mfields_zero(flds, JYI);
  psc_mfields_zero(flds, JZI);

  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = ppsc->patch + p;
    fields_t *pf = psc_mfields_get_patch(flds, p);
    fields_single_t fld;
    // FIXME, can do -1 .. 1?
    int ib[3] = { 0, -2, -2 };
    int ie[3] = { 1, patch->ldims[1] + 2, patch->ldims[2] + 2 };
    fields_single_alloc(&fld, ib, ie, 9, 0); // JXI .. HZ
    for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
      for (int iy = -2; iy < patch->ldims[1] + 2; iy++) {
	F3_S(&fld, EX, 0,iy,iz) = F3(pf, EX, 0,iy,iz);
	F3_S(&fld, EY, 0,iy,iz) = F3(pf, EY, 0,iy,iz);
	F3_S(&fld, EZ, 0,iy,iz) = F3(pf, EZ, 0,iy,iz);
	F3_S(&fld, HX, 0,iy,iz) = F3(pf, HX, 0,iy,iz);
	F3_S(&fld, HY, 0,iy,iz) = F3(pf, HY, 0,iy,iz);
	F3_S(&fld, HZ, 0,iy,iz) = F3(pf, HZ, 0,iy,iz);
      }
    }

    do_push_part_1vb_yz(p, &fld, psc_mparticles_get_patch(particles, p));

    for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
      for (int iy = -2; iy < patch->ldims[1] + 2; iy++) {
	F3(pf, JXI, 0,iy,iz) += F3_S(&fld, JXI, 0,iy,iz);
	F3(pf, JYI, 0,iy,iz) += F3_S(&fld, JYI, 0,iy,iz);
	F3(pf, JZI, 0,iy,iz) += F3_S(&fld, JZI, 0,iy,iz);
      }
    }
    fields_single_free(&fld);
  }
  prof_stop(pr);

  psc_mfields_put_cf(flds, flds_base, JXI, JXI + 3);
  psc_mparticles_put_cf(particles, particles_base);
}

