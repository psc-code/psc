
#include "psc_output_fields_item_private.h"
#include <psc_bnd.h>
#include <fields.hxx>

#include <math.h>

using Fields = Fields3d<fields_t>;
using real_t = mparticles_t::real_t;

#include "common_moments.cxx"

// ======================================================================
// boundary stuff FIXME, should go elsewhere...

static void
add_ghosts_reflecting_lo(fields_t flds, int p, int d, int mb, int me)
{
  Fields F(flds);
  struct psc_patch *patch = ppsc->patch + p;

  int bx = patch->ldims[0] == 1 ? 0 : 1;
  if (d == 1) {
    for (int iz = -1; iz < patch->ldims[2] + 1; iz++) {
      for (int ix = -bx; ix < patch->ldims[0] + bx; ix++) {
	int iy = 0; {
	  for (int m = mb; m < me; m++) {
	    F(m, ix,iy,iz) += F(m, ix,iy-1,iz);
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = 0*-1; iy < patch->ldims[1] + 0*1; iy++) {
      for (int ix = -bx; ix < patch->ldims[0] + bx; ix++) {
	int iz = 0; {
	  for (int m = mb; m < me; m++) {
	    F(m, ix,iy,iz) += F(m, ix,iy,iz-1);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_reflecting_hi(fields_t flds, int p, int d, int mb, int me)
{
  Fields F(flds);
  struct psc_patch *patch = ppsc->patch + p;

  int bx = patch->ldims[0] == 1 ? 0 : 1;
  if (d == 1) {
    for (int iz = -1; iz < patch->ldims[2] + 1; iz++) {
      for (int ix = -bx; ix < patch->ldims[0] + bx; ix++) {
	int iy = patch->ldims[1] - 1; {
	  for (int m = mb; m < me; m++) {
	    F(m, ix,iy,iz) += F(m, ix,iy+1,iz);
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = 0*-1; iy < patch->ldims[1] + 0*1; iy++) {
      for (int ix = -bx; ix < patch->ldims[0] + bx; ix++) {
	int iz = patch->ldims[2] - 1; {
	  for (int m = mb; m < me; m++) {
	    F(m, ix,iy,iz) += F(m, ix,iy,iz+1);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_boundary(fields_t res, int p, int mb, int me)
{
  // lo
  for (int d = 0; d < 3; d++) {
    if (psc_at_boundary_lo(ppsc, p, d)) {
      if (ppsc->domain.bnd_part_lo[d] == BND_PART_REFLECTING ||
	  ppsc->domain.bnd_part_lo[d] == BND_PART_OPEN) {
	add_ghosts_reflecting_lo(res, p, d, mb, me);
      }
    }
  }
  // hi
  for (int d = 0; d < 3; d++) {
    if (psc_at_boundary_hi(ppsc, p, d)) {
      if (ppsc->domain.bnd_part_hi[d] == BND_PART_REFLECTING ||
	  ppsc->domain.bnd_part_hi[d] == BND_PART_OPEN) {
	add_ghosts_reflecting_hi(res, p, d, mb, me);
      }
    }
  }
}

static void
run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	struct psc_mparticles *mprts_base, struct psc_mfields *mres,
	void (*do_run)(int p, fields_t res, mparticles_t& mprts))
{
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  mfields_t mf_res(mres);
  
  for (int p = 0; p < mprts->n_patches(); p++) {
    mf_res[p].zero();
    do_run(p, mf_res[p], mprts);
    add_ghosts_boundary(mf_res[p], p, 0, mres->nr_fields);
  }

  mprts.put_as(mprts_base, MP_DONT_COPY);
}

// ======================================================================
// n_1st

static void
do_n_1st_run(int p, fields_t flds, mparticles_t& mprts)
{
  const Grid_t& grid = ppsc->grid();
  mparticles_t::patch_t& prts = mprts[p];
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int m = prt->kind();
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, m, 1.f);
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
do_v_1st_run(int p, fields_t flds, mparticles_t& mprts)
{
  const Grid_t& grid = ppsc->grid();
  mparticles_t::patch_t& prts = mprts[p];
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 3;

    real_t vxi[3];
    particle_calc_vxi(prt, vxi);

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, vxi[m]);
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
do_p_1st_run(int p, fields_t flds, mparticles_t& mprts)
{
  const Grid_t& grid = ppsc->grid();
  mparticles_t::patch_t& prts = mprts[p];
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 3;
    real_t *pxi = &prt->pxi;

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, mprts->prt_mni(*prt) * pxi[m]);
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
do_vv_1st_run(int p, fields_t flds, mparticles_t& mprts)
{
  const Grid_t& grid = ppsc->grid();
  mparticles_t::patch_t& prts = mprts[p];
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 3;

    real_t vxi[3];
    particle_calc_vxi(prt, vxi);

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, vxi[m] * vxi[m]);
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
do_T_1st_run(int p, fields_t flds, mparticles_t& mprts)
{
  const Grid_t& grid = ppsc->grid();
  mparticles_t::patch_t& prts = mprts[p];
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 6;

    real_t vxi[3];
    particle_calc_vxi(prt, vxi);
    real_t *pxi = &prt->pxi;
    real_t vx[3];
    vx[0] = vxi[0] * cos(ppsc->prm.theta_xz) - vxi[2] * sin(ppsc->prm.theta_xz);
    vx[1] = vxi[1];
    vx[2] = vxi[0] * sin(ppsc->prm.theta_xz) + vxi[2] * cos(ppsc->prm.theta_xz);
    real_t px[3];
    px[0] = pxi[0] * cos(ppsc->prm.theta_xz) - pxi[2] * sin(ppsc->prm.theta_xz);
    px[1] = pxi[1];
    px[2] = pxi[0] * sin(ppsc->prm.theta_xz) + pxi[2] * cos(ppsc->prm.theta_xz);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 0, mprts->prt_mni(*prt) * px[0] * vx[0]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 1, mprts->prt_mni(*prt) * px[1] * vx[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 2, mprts->prt_mni(*prt) * px[2] * vx[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 3, mprts->prt_mni(*prt) * px[0] * vx[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4, mprts->prt_mni(*prt) * px[0] * vx[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 5, mprts->prt_mni(*prt) * px[1] * vx[2]);
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
do_Tvv_1st_run(int p, fields_t flds, mparticles_t& mprts)
{
  const Grid_t& grid = ppsc->grid();
  mparticles_t::patch_t& prts = mprts[p];
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 6;

    real_t vxi[3];
    particle_calc_vxi(prt, vxi);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 0, mprts->prt_mni(*prt) * vxi[0] * vxi[0]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 1, mprts->prt_mni(*prt) * vxi[1] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 2, mprts->prt_mni(*prt) * vxi[2] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 3, mprts->prt_mni(*prt) * vxi[0] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4, mprts->prt_mni(*prt) * vxi[0] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 5, mprts->prt_mni(*prt) * vxi[1] * vxi[2]);
  }
}

static void
Tvv_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds_base, mprts_base, mres, do_Tvv_1st_run);
}

#if 0

// ======================================================================
// nvt_1st

static void
do_nvt_a_1st_run(int p, fields_t flds, mparticles_t::patch_t& prts)
{
  const Grid_t& grid = ppsc->grid();
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 10;

    real_t vxi[3];
    particle_calc_vxi(prt, vxi);

    // density
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm, 1.f);
    // velocity
    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m + 1, vxi[m]);
    }
  }
}

static void
do_nvt_b_1st_run(int p, fields_t flds, mparticles_t::patch_t& prts)
{
  Fields F(flds);
  const Grid_t& grid = ppsc->grid();
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 10;

    real_t *xi = &prt->xi;
    real_t u = xi[0] * dxi - .5;
    real_t v = xi[1] * dyi - .5;
    real_t w = xi[2] * dzi - .5;
    int jx = fint(u);
    int jy = fint(v);
    int jz = fint(w);
    real_t h1 = u - jx;
    real_t h2 = v - jy;
    real_t h3 = w - jz;

    real_t g0x = 1.f - h1;
    real_t g0y = 1.f - h2;
    real_t g0z = 1.f - h3;
    real_t g1x = h1;
    real_t g1y = h2;
    real_t g1z = h3;

    int jxd = 1, jyd = 1, jzd = 1;
    if (ppsc->domain.gdims[0] == 1) {
      jx = 0; g0x = 1.; g1x = 0.; jxd = 0;
    }
    if (ppsc->domain.gdims[1] == 1) {
      jy = 0; g0y = 1.; g1y = 0.; jyd = 0;
    }
    if (ppsc->domain.gdims[2] == 1) {
      jz = 0; g0z = 1.; g1z = 0.; jzd = 0;
    }

    assert(jx >= -1 && jx < patch->ldims[0]);
    assert(jy >= -1 && jy < patch->ldims[1]);
    assert(jz >= -1 && jz < patch->ldims[2]);

    real_t vxi[3];
    particle_calc_vxi(prt, vxi);
    for (int d = 0; d < 3; d++) {
      int m = mm + 1 + d;
      double vavg = (g0x*g0y*g0z * F(m, jx    ,jy    ,jz    ) +
		     g1x*g0y*g0z * F(m, jx+jxd,jy    ,jz    ) +
		     g0x*g1y*g0z * F(m, jx    ,jy+jyd,jz    ) +
		     g1x*g1y*g0z * F(m, jx+jxd,jy+jyd,jz    ) +
		     g0x*g0y*g1z * F(m, jx    ,jy    ,jz+jzd) +
		     g1x*g0y*g1z * F(m, jx+jxd,jy    ,jz+jzd) +
		     g0x*g1y*g1z * F(m, jx    ,jy+jyd,jz+jzd) +
		     g1x*g1y*g1z * F(m, jx+jxd,jy+jyd,jz+jzd));
      vxi[d] -= vavg;
    }
    real_t *pxi = vxi;
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 0, pxi[0] * vxi[0]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 1, pxi[1] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 2, pxi[2] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 3, pxi[0] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 4, pxi[0] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 5, pxi[1] * vxi[2]);
  }
}

static void
do_nvp_1st_run(int p, fields_t flds, mparticles_t::patch_t& prts)
{
  Fields F(flds);
  const Grid_t& grid = ppsc->grid();
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 10;

    real_t *xi = &prt->xi;
    real_t u = xi[0] * dxi - .5;
    real_t v = xi[1] * dyi - .5;
    real_t w = xi[2] * dzi - .5;
    int jx = fint(u);
    int jy = fint(v);
    int jz = fint(w);

    assert(jx >= -1 && jx < patch->ldims[0]);
    assert(jy >= -1 && jy < patch->ldims[1]);
    assert(jz >= -1 && jz < patch->ldims[2]);

    real_t vxi[3];
    particle_calc_vxi(prt, vxi);
    real_ty *pxi = vxi;
    // density
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm, 1.f);
    // velocity
    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m + 1, vxi[m]);
    }
    // pi * vj
    // FIXME, this is really vxi * vxi, etc
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 0, pxi[0] * vxi[0]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 1, pxi[1] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 2, pxi[2] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 3, pxi[0] * vxi[1]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 4, pxi[1] * vxi[2]);
    DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4 + 5, pxi[2] * vxi[0]);
  }
}

static void
nvt_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
		struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  const Grid_t& grid = ppsc->grid();
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  mfields_t mf_res(mres);

  for (int p = 0; p < mres->nr_patches; p++) {
    mf_res[p].zero();
    do_nvt_a_1st_run(p, mf_res[p], mprts[p].range());
    add_ghosts_boundary(mf_res[p], p, 0, mres->nr_fields);
  }

  psc_bnd_add_ghosts(item->bnd, mres, 0, mres->nr_fields);
  psc_bnd_fill_ghosts(item->bnd, mres, 0, mres->nr_fields);

  for (int p = 0; p < mres->nr_patches; p++) {
    Fields R(mf_res[p]);

    // fix up zero density cells
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      psc_foreach_3d(ppsc, p, ix, iy, iz, 1, 1) {
	if (R(10*m, ix,iy,iz) == 0.0) {
	  R(10*m, ix,iy,iz) = 0.00001;
	} psc_foreach_3d_end;
      }
    }    

    // normalize v moments
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      for (int mm = 0; mm < 3; mm++) {
	psc_foreach_3d(ppsc, p, ix, iy, iz, 1, 1) {
	  R(10*m + mm + 1, ix,iy,iz) /= R(10*m, ix,iy,iz);
	} psc_foreach_3d_end;
      }
    }

    // calculate <(v-U)(v-U)> moments
    do_nvt_b_1st_run(p, mf_res[p], mprts[p].range());
  }

  mprts.put_as(mprts_base, MP_DONT_COPY);

  for (int m = 0; m < ppsc->nr_kinds; m++) {
    psc_bnd_add_ghosts(item->bnd, mres, 10*m + 4, 10*m + 10);
  }
  psc_bnd_fill_ghosts(item->bnd, mres, 0, mres->nr_fields);
}

static void
nvp_1st_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
		struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  mfields_t mf_res(mres);

  for (int p = 0; p < mres->nr_patches; p++) {
    mf_res[p].zero();
    do_nvp_1st_run(p, mf_res[p], mprts[p].range());
    add_ghosts_boundary(mf_res[p], p, 0, mres->nr_fields);
  }

  psc_bnd_add_ghosts(item->bnd, mres, 0, mres->nr_fields);
  psc_bnd_fill_ghosts(item->bnd, mres, 0, mres->nr_fields);

  const int mm2mx[6] = { 0, 1, 2, 0, 1, 2 };
  const int mm2my[6] = { 0, 1, 2, 1, 2, 0 };

  for (int p = 0; p < mres->nr_patches; p++) {
    Fields R(mf_res[p]);

    // fix up zero density cells
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      foreach_3d(ppsc, p, ix, iy, iz, 1, 1) {
	if (R(10*m, ix,iy,iz) == 0.0) {
	  R(10*m, ix,iy,iz) = 0.00001;
	} foreach_3d_end;
      }
    }    

    // normalize v moments
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      for (int mm = 0; mm < 3; mm++) {
	foreach_3d(ppsc, p, ix, iy, iz, 1, 1) {
	  R(10*m + mm + 1, ix,iy,iz) /= R(10*m, ix,iy,iz);
	} foreach_3d_end;
      }
    }

    // calculate <(v - U)(v - U)> moments
    for (int m = 0; m < ppsc->nr_kinds; m++) {
      for (int mm = 0; mm < 6; mm++) {
	int mx = mm2mx[mm], my = mm2my[mm];
	foreach_3d(ppsc, p, ix, iy, iz, 1, 1) {
	  R(10*m + 4 + mm, ix,iy,iz) =
	    R(10*m + 4 + mm, ix,iy,iz) / R(10*m, ix,iy,iz) - 
	    R(10*m + 1 + mx, ix,iy,iz) * R(10*m + 1 + my, ix,iy,iz);
	} foreach_3d_end;
      }
    }
  }

  mprts.put_as(mprts_base, MP_DONT_COPY);
}

#endif

// ======================================================================

#define MAKE_OP1(TYPE, NAME, FNAME, RUN)				\
  struct psc_output_fields_item_ops_##NAME##TYPE : psc_output_fields_item_ops { \
    psc_output_fields_item_ops_##NAME##TYPE() {				\
    name               = #NAME #TYPE;					\
    nr_comp	       = 1;						\
    fld_names[0]       = FNAME;						\
    run_all            = RUN;						\
    flags              = POFI_ADD_GHOSTS | POFI_BY_KIND;		\
  }									\
  } psc_output_fields_item_##NAME##TYPE##_ops;

#define MAKE_OP3(TYPE, NAME, FNAMEX, FNAMEY, FNAMEZ, RUN)		\
  struct psc_output_fields_item_ops_##NAME##TYPE : psc_output_fields_item_ops { \
    psc_output_fields_item_ops_##NAME##TYPE() {				\
    name               = #NAME #TYPE;					\
    nr_comp	       = 3;						\
    fld_names[0]       = FNAMEX;					\
    fld_names[1]       = FNAMEY;					\
    fld_names[2]       = FNAMEZ;					\
    run_all            = RUN;						\
    flags              = POFI_ADD_GHOSTS | POFI_BY_KIND;		\
  }									\
  } psc_output_fields_item_##NAME##TYPE##_ops;

#define MAKE_OP6(TYPE, NAME, FNAMEX, FNAMEY, FNAMEZ, FNAME3, FNAME4, FNAME5, RUN) \
  struct psc_output_fields_item_ops_##NAME##TYPE : psc_output_fields_item_ops { \
    psc_output_fields_item_ops_##NAME##TYPE() {				\
    name               = #NAME #TYPE;					\
    nr_comp	       = 6;						\
    fld_names[0]       = FNAMEX;					\
    fld_names[1]       = FNAMEY;					\
    fld_names[2]       = FNAMEZ;					\
    fld_names[3]       = FNAME3;					\
    fld_names[4]       = FNAME4;					\
    fld_names[5]       = FNAME5;					\
    run_all            = RUN;						\
    flags              = POFI_ADD_GHOSTS | POFI_BY_KIND;		\
  }									\
  } psc_output_fields_item_##NAME##TYPE##_ops;


#define MAKE_POFI_OPS(TYPE)						\
  MAKE_OP1(TYPE, n_1st_, "n", n_1st_run_all)				\
  MAKE_OP3(TYPE, v_1st_, "vx", "vy", "vz", v_1st_run_all)		\
  MAKE_OP3(TYPE, p_1st_, "px", "py", "pz", p_1st_run_all)		\
  MAKE_OP3(TYPE, vv_1st_, "vxvx", "vyvy", "vzvz", vv_1st_run_all)	\
  MAKE_OP6(TYPE, T_1st_, "Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz", T_1st_run_all)	\
  MAKE_OP6(TYPE, Tvv_1st_, "vxvx", "vyvy", "vzvz", "vxvy", "vxvz", "vyvz", Tvv_1st_run_all)	\

#if 0
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

#endif

