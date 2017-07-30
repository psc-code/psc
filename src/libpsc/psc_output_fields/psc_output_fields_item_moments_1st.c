
#include "psc_output_fields_item_private.h"
#include <psc_bnd.h>

#include <math.h>

#include "common_moments.c"

// ======================================================================
// boundary stuff FIXME, should go elsewhere...

static void
add_ghosts_reflecting_lo(struct psc_fields *pf, int d, int mb, int me)
{
  struct psc_patch *patch = ppsc->patch + pf->p;

  int bx = patch->ldims[0] == 1 ? 0 : 1;
  if (d == 1) {
    for (int iz = -1; iz < patch->ldims[2] + 1; iz++) {
      for (int ix = -bx; ix < patch->ldims[0] + bx; ix++) {
	int iy = 0; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy-1,iz);
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = -1; iy < patch->ldims[1] + 1; iy++) {
      for (int ix = -bx; ix < patch->ldims[0] + bx; ix++) {
	int iz = 0; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy,iz-1);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_reflecting_hi(struct psc_fields *pf, int d, int mb, int me)
{
  struct psc_patch *patch = ppsc->patch + pf->p;

  int bx = patch->ldims[0] == 1 ? 0 : 1;
  if (d == 1) {
    for (int iz = -1; iz < patch->ldims[2] + 1; iz++) {
      for (int ix = -bx; ix < patch->ldims[0] + bx; ix++) {
	int iy = patch->ldims[1] - 1; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy+1,iz);
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = -1; iy < patch->ldims[1] + 1; iy++) {
      for (int ix = -bx; ix < patch->ldims[0] + bx; ix++) {
	int iz = patch->ldims[2] - 1; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy,iz+1);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_boundary(struct psc_fields *res, int mb, int me)
{
  // lo
  for (int d = 0; d < 3; d++) {
    if (ppsc->patch[res->p].off[d] == 0) {
      if (ppsc->domain.bnd_part_lo[d] == BND_PART_REFLECTING ||
	  ppsc->domain.bnd_part_lo[d] == BND_PART_OPEN) {
	add_ghosts_reflecting_lo(res, d, mb, me);
      }
    }
  }
  // hi
  for (int d = 0; d < 3; d++) {
    if (ppsc->patch[res->p].off[d] + ppsc->patch[res->p].ldims[d] == ppsc->domain.gdims[d]) {
      if (ppsc->domain.bnd_part_hi[d] == BND_PART_REFLECTING ||
	  ppsc->domain.bnd_part_hi[d] == BND_PART_OPEN) {
	add_ghosts_reflecting_hi(res, d, mb, me);
      }
    }
  }
}

static void
run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	struct psc_mparticles *mprts_base, struct psc_mfields *mres,
	void (*do_run)(int p, struct psc_fields *flds,
		       particle_iter_t prt_begin, particle_iter_t prt_end))
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_fields *res = psc_mfields_get_patch(mres, p);
    psc_particles_reorder(prts); // FIXME
    psc_fields_zero_range(res, 0, res->nr_comp);
    do_run(res->p, res, particle_iter_begin_prts(prts), particle_iter_end_prts(prts));
    add_ghosts_boundary(res, 0, res->nr_comp);
  }

  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
}

// ======================================================================
// n_1st

static void
do_n_1st_run(int p, fields_t *pf, particle_iter_t prt_begin, particle_iter_t prt_end)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];


  PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int m = particle_kind(prt);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, m, 1.f);
  }
}

static void
n_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	      struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds_base, mprts_base, mres, do_n_1st_run);
}

// ======================================================================
// v_1st

static void
do_v_1st_run(int p, fields_t *pf, particle_iter_t prt_begin, particle_iter_t prt_end)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int mm = particle_kind(prt) * 3;

    particle_real_t vxi[3];
    particle_calc_vxi(prt, vxi);

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + m, vxi[m]);
    }
  }
}

static void
v_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	      struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds_base, mprts_base, mres, do_v_1st_run);
}

// ======================================================================
// p_1st

static void
do_p_1st_run(int p, fields_t *pf, particle_iter_t prt_begin, particle_iter_t prt_end)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int mm = particle_kind(prt) * 3;
    particle_real_t *pxi = &prt->pxi;

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + m, particle_mni(prt) * pxi[m]);
    }
  }
}

static void
p_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	      struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds_base, mprts_base, mres, do_p_1st_run);
}

// ======================================================================
// vv_1st

static void
do_vv_1st_run(int p, fields_t *pf, particle_iter_t prt_begin, particle_iter_t prt_end)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int mm = particle_kind(prt) * 3;

    particle_real_t vxi[3];
    particle_calc_vxi(prt, vxi);

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + m, vxi[m] * vxi[m]);
    }
  }
}

static void
vv_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	       struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds_base, mprts_base, mres, do_vv_1st_run);
}

// ======================================================================
// T_1st

static void
do_T_1st_run(int p, fields_t *pf, particle_iter_t prt_begin, particle_iter_t prt_end)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int mm = particle_kind(prt) * 6;

    particle_real_t vxi[3];
    particle_calc_vxi(prt, vxi);
    particle_real_t *pxi = &prt->pxi;
    particle_real_t vx[3] = {
      vxi[0] * cos(ppsc->prm.theta_xz) - vxi[2] * sin(ppsc->prm.theta_xz),
      vxi[1],
      vxi[0] * sin(ppsc->prm.theta_xz) + vxi[2] * cos(ppsc->prm.theta_xz),
    };
    particle_real_t px[3] = {
      pxi[0] * cos(ppsc->prm.theta_xz) - pxi[2] * sin(ppsc->prm.theta_xz),
      pxi[1],
      pxi[0] * sin(ppsc->prm.theta_xz) + pxi[2] * cos(ppsc->prm.theta_xz),
    };
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 0, particle_mni(prt) * px[0] * vx[0]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 1, particle_mni(prt) * px[1] * vx[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 2, particle_mni(prt) * px[2] * vx[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 3, particle_mni(prt) * px[0] * vx[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 4, particle_mni(prt) * px[0] * vx[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 5, particle_mni(prt) * px[1] * vx[2]);
  }
}

static void
T_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	      struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds_base, mprts_base, mres, do_T_1st_run);
}

// ======================================================================
// Tvv_1st

static void
do_Tvv_1st_run(int p, fields_t *pf, particle_iter_t prt_begin, particle_iter_t prt_end)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int mm = particle_kind(prt) * 6;

    particle_real_t vxi[3];
    particle_calc_vxi(prt, vxi);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 0, particle_mni(prt) * vxi[0] * vxi[0]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 1, particle_mni(prt) * vxi[1] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 2, particle_mni(prt) * vxi[2] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 3, particle_mni(prt) * vxi[0] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 4, particle_mni(prt) * vxi[0] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 5, particle_mni(prt) * vxi[1] * vxi[2]);
  }
}

static void
Tvv_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds_base, mprts_base, mres, do_Tvv_1st_run);
}

// ======================================================================
// nvt_1st

static void
do_nvt_a_1st_run(int p, fields_t *pf, particle_iter_t prt_begin, particle_iter_t prt_end)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int mm = particle_kind(prt) * 10;

    particle_real_t vxi[3];
    particle_calc_vxi(prt, vxi);

    // density
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm, 1.f);
    // velocity
    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + m + 1, vxi[m]);
    }
  }
}

static void
do_nvt_b_1st_run(int p, fields_t *pf, particle_iter_t prt_begin, particle_iter_t prt_end)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end) {
    particle_t *prt = particle_iter_deref(prt_iter);
    int mm = particle_kind(prt) * 10;

    particle_real_t *xi = &prt->xi;					\
    particle_real_t u = xi[0] * dxi - .5;				\
    particle_real_t v = xi[1] * dyi - .5;				\
    particle_real_t w = xi[2] * dzi - .5;				\
    int jx = particle_real_fint(u);					\
    int jy = particle_real_fint(v);					\
    int jz = particle_real_fint(w);					\
    particle_real_t h1 = u - jx;					\
    particle_real_t h2 = v - jy;					\
    particle_real_t h3 = w - jz;					\
    									\
    particle_real_t g0x = 1.f - h1;					\
    particle_real_t g0y = 1.f - h2;					\
    particle_real_t g0z = 1.f - h3;					\
    particle_real_t g1x = h1;						\
    particle_real_t g1y = h2;						\
    particle_real_t g1z = h3;						\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (ppsc->domain.gdims[0] == 1) {					\
      jx = 0; g0x = 1.; g1x = 0.; jxd = 0;				\
    }									\
    if (ppsc->domain.gdims[1] == 1) {					\
      jy = 0; g0y = 1.; g1y = 0.; jyd = 0;				\
    }									\
    if (ppsc->domain.gdims[2] == 1) {					\
      jz = 0; g0z = 1.; g1z = 0.; jzd = 0;				\
    }									\
    									\
    assert(jx >= -1 && jx < patch->ldims[0]);				\
    assert(jy >= -1 && jy < patch->ldims[1]);				\
    assert(jz >= -1 && jz < patch->ldims[2]);				\
    									\
    particle_real_t vxi[3];
    particle_calc_vxi(prt, vxi);
    for (int d = 0; d < 3; d++) {
      int m = mm + 1 + d;
      double vavg = (g0x*g0y*g0z * F3(pf, m, jx    ,jy    ,jz    ) +
		     g1x*g0y*g0z * F3(pf, m, jx+jxd,jy    ,jz    ) +
		     g0x*g1y*g0z * F3(pf, m, jx    ,jy+jyd,jz    ) +
		     g1x*g1y*g0z * F3(pf, m, jx+jxd,jy+jyd,jz    ) +
		     g0x*g0y*g1z * F3(pf, m, jx    ,jy    ,jz+jzd) +
		     g1x*g0y*g1z * F3(pf, m, jx+jxd,jy    ,jz+jzd) +
		     g0x*g1y*g1z * F3(pf, m, jx    ,jy+jyd,jz+jzd) +
		     g1x*g1y*g1z * F3(pf, m, jx+jxd,jy+jyd,jz+jzd));
      vxi[d] -= vavg;
    }
    particle_real_t *pxi = vxi;
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 4 + 0, pxi[0] * vxi[0]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 4 + 1, pxi[1] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 4 + 2, pxi[2] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 4 + 3, pxi[0] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 4 + 4, pxi[0] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, pf, mm + 4 + 5, pxi[1] * vxi[2]);
  }
}

static void
do_nvp_1st_run(int p, fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int mm = particle_kind(part) * 10;

    particle_real_t *xi = &part->xi;					\
    particle_real_t u = xi[0] * dxi - .5;				\
    particle_real_t v = xi[1] * dyi - .5;				\
    particle_real_t w = xi[2] * dzi - .5;				\
    int jx = particle_real_fint(u);					\
    int jy = particle_real_fint(v);					\
    int jz = particle_real_fint(w);					\
    									\
    assert(jx >= -1 && jx < patch->ldims[0]);				\
    assert(jy >= -1 && jy < patch->ldims[1]);				\
    assert(jz >= -1 && jz < patch->ldims[2]);				\
    									\
    particle_real_t vxi[3];
    particle_calc_vxi(part, vxi);
    particle_real_t *pxi = vxi;
    // density
    DEPOSIT_TO_GRID_1ST_CC(part, pf, mm, 1.f);
    // velocity
    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + m + 1, vxi[m]);
    }
    // pi * vj
    // FIXME, this is really vxi * vxi, etc
    DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + 4 + 0, pxi[0] * vxi[0]);
    DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + 4 + 1, pxi[1] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + 4 + 2, pxi[2] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + 4 + 3, pxi[0] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + 4 + 4, pxi[1] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(part, pf, mm + 4 + 5, pxi[2] * vxi[0]);
  }
}

static void
nvt_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
		struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);

  for (int p = 0; p < mres->nr_patches; p++) {
    struct psc_fields *res = psc_mfields_get_patch(mres, p);
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    psc_fields_zero_range(res, 0, res->nr_comp);
    do_nvt_a_1st_run(res->p, res, particle_iter_begin_prts(prts),
		     particle_iter_end_prts(prts));
    add_ghosts_boundary(res, 0, res->nr_comp);
  }

  psc_bnd_add_ghosts(item->bnd, mres, 0, mres->nr_fields);
  psc_bnd_fill_ghosts(item->bnd, mres, 0, mres->nr_fields);

  for (int p = 0; p < mres->nr_patches; p++) {
    struct psc_fields *res = psc_mfields_get_patch(mres, p);
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);

    // fix up zero density cells
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      psc_foreach_3d(ppsc, res->p, ix, iy, iz, 1, 1) {
	if (F3(res, 10*m, ix,iy,iz) == 0.0) {
	  F3(res, 10*m, ix,iy,iz) = 0.00001;
	} psc_foreach_3d_end;
      }
    }    

    // normalize v moments
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      for (int mm = 0; mm < 3; mm++) {
	psc_foreach_3d(ppsc, res->p, ix, iy, iz, 1, 1) {
	  F3(res, 10*m + mm + 1, ix,iy,iz) /= F3(res, 10*m, ix,iy,iz);
	} psc_foreach_3d_end;
      }
    }

    // calculate <(v-U)(v-U)> moments
    do_nvt_b_1st_run(res->p, res, particle_iter_begin_prts(prts),
		     particle_iter_end_prts(prts));
  }

  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);

  for (int m = 0; m < ppsc->nr_kinds; m++) {
    psc_bnd_add_ghosts(item->bnd, mres, 10*m + 4, 10*m + 10);
  }
  psc_bnd_fill_ghosts(item->bnd, mres, 0, mres->nr_fields);
}

static void
nvp_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
		struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);

  for (int p = 0; p < mres->nr_patches; p++) {
    struct psc_fields *res = psc_mfields_get_patch(mres, p);
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    psc_fields_zero_range(res, 0, res->nr_comp);
    do_nvp_1st_run(res->p, res, prts);
    add_ghosts_boundary(res, 0, res->nr_comp);
  }

  psc_bnd_add_ghosts(item->bnd, mres, 0, mres->nr_fields);
  psc_bnd_fill_ghosts(item->bnd, mres, 0, mres->nr_fields);

  const int mm2mx[6] = { 0, 1, 2, 0, 1, 2 };
  const int mm2my[6] = { 0, 1, 2, 1, 2, 0 };

  for (int p = 0; p < mres->nr_patches; p++) {
    struct psc_fields *res = psc_mfields_get_patch(mres, p);

    // fix up zero density cells
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      foreach_3d(ppsc, res->p, ix, iy, iz, 1, 1) {
	if (F3(res, 10*m, ix,iy,iz) == 0.0) {
	  F3(res, 10*m, ix,iy,iz) = 0.00001;
	} foreach_3d_end;
      }
    }    

    // normalize v moments
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      for (int mm = 0; mm < 3; mm++) {
	foreach_3d(ppsc, res->p, ix, iy, iz, 1, 1) {
	  F3(res, 10*m + mm + 1, ix,iy,iz) /= F3(res, 10*m, ix,iy,iz);
	} foreach_3d_end;
      }
    }

    // calculate <(v - U)(v - U)> moments
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      for (int mm = 0; mm < 6; mm++) {
	int mx = mm2mx[mm], my = mm2my[mm];
	foreach_3d(ppsc, res->p, ix, iy, iz, 1, 1) {
	  F3(res, 10*m + 4 + mm, ix,iy,iz) =
	    F3(res, 10*m + 4 + mm, ix,iy,iz) / F3(res, 10*m, ix,iy,iz) - 
	    F3(res, 10*m + 1 + mx, ix,iy,iz) * F3(res, 10*m + 1 + my, ix,iy,iz);
	} foreach_3d_end;
      }
    }
  }

  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
}

// ======================================================================

#define MAKE_POFI_OPS(TYPE)						\
struct psc_output_fields_item_ops psc_output_fields_item_n_1st_##TYPE##_ops = { \
  .name               = "n_1st_" #TYPE,					\
  .nr_comp	      = 1,						\
  .fld_names	      = { "n" },					\
  .run_all            = n_1st_run_all,					\
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,			\
};									\
									\
struct psc_output_fields_item_ops psc_output_fields_item_v_1st_##TYPE##_ops = { \
  .name               = "v_1st_" #TYPE,					\
  .nr_comp	      = 3,						\
  .fld_names	      = { "vx", "vy", "vz" },				\
  .run_all            = v_1st_run_all,					\
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,			\
};									\
									\
struct psc_output_fields_item_ops psc_output_fields_item_p_1st_##TYPE##_ops = { \
  .name               = "p_1st_" #TYPE,					\
  .nr_comp	      = 3,						\
  .fld_names	      = { "px", "py", "pz" },				\
  .run_all            = p_1st_run_all,					\
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,			\
};									\
									\
struct psc_output_fields_item_ops psc_output_fields_item_vv_1st_##TYPE##_ops = { \
  .name               = "vv_1st_" #TYPE,				\
  .nr_comp	      = 3,						\
  .fld_names	      = { "vxvx", "vyvy", "vzvz" },			\
  .run_all            = vv_1st_run_all,					\
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,			\
};									\
									\
struct psc_output_fields_item_ops psc_output_fields_item_T_1st_##TYPE##_ops = { \
  .name               = "T_1st_" #TYPE,					\
  .nr_comp	      = 6,						\
  .fld_names	      = { "Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz" },	\
  .run_all            = T_1st_run_all,					\
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,			\
};									\
									\
struct psc_output_fields_item_ops psc_output_fields_item_Tvv_1st_##TYPE##_ops = { \
  .name               = "Tvv_1st_" #TYPE,					\
  .nr_comp	      = 6,						\
  .fld_names	      = { "vxvx", "vyvy", "vzvz", "vxvy", "vxvz", "vyvz" },	\
  .run_all            = Tvv_1st_run_all,					\
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,			\
};									\
									\
struct psc_output_fields_item_ops psc_output_fields_item_nvt_1st_##TYPE##_ops = { \
  .name               = "nvt_1st_" #TYPE,				\
  .nr_comp	      = 10,						\
  .fld_names	      = { "n", "vx", "vy", "vz",			\
			  "Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz" },	\
  .run_all            = nvt_1st_run_all,				\
  .flags              = POFI_BY_KIND,					\
};									\
									\
 /* FIXME, should include mass in Tij?! */				\
struct psc_output_fields_item_ops psc_output_fields_item_nvp_1st_##TYPE##_ops = { \
  .name               = "nvp_1st_" #TYPE,				\
  .nr_comp	      = 10,						\
  .fld_names	      = { "n", "vx", "vy", "vz",			\
			  "Txx", "Tyy", "Tzz", "Txy", "Tyz", "Tzx" },	\
  .run_all            = nvp_1st_run_all,				\
  .flags              = POFI_BY_KIND,					\
};									\
									\

