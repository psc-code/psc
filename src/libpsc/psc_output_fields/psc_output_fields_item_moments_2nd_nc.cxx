
#include "psc_output_fields_item_private.h"
#include "fields.hxx"

using Fields = Fields3d<fields_t>;

#include <math.h>

#include "common_moments.cxx"

// ======================================================================
// boundary stuff FIXME, should go elsewhere, and FIXME, duplicated, kinda,
// from cell-centered 1st order

static void
add_ghosts_reflecting_lo(int p, fields_t flds, int d, int mb, int me)
{
  Fields F(flds);
  struct psc_patch *patch = ppsc->patch + p;

  if (d == 1) {
    for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iy = 0; {
	  for (int m = mb; m < me; m++) {
	    F(m, ix,iy+1,iz) += F(m, ix,iy-1,iz);
	    F(m, ix,iy-1,iz) = 0.;
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = 0; iy < patch->ldims[1]; iy++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iz = 0; {
	  for (int m = mb; m < me; m++) {
	    F(m, ix,iy,iz+1) += F(m, ix,iy,iz-1);
	    F(m, ix,iy,iz-1) = 0.;
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_reflecting_hi(int p, fields_t flds, int d, int mb, int me)
{
  Fields F(flds);
  struct psc_patch *patch = ppsc->patch + p;

  if (d == 1) {
    for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iy = patch->ldims[1]; {
	  for (int m = mb; m < me; m++) {
	    F(m, ix,iy-1,iz) += F(m, ix,iy+1,iz);
	    F(m, ix,iy+1,iz) = 0.;
	  }
	}
      }
    }
  } else if (d == 2) {
    for (int iy = 0; iy < patch->ldims[1]; iy++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iz = patch->ldims[2]; {
	  for (int m = mb; m < me; m++) {
	    F(m, ix,iy,iz-1) += F(m, ix,iy,iz+1);
	    F(m, ix,iy,iz+1) += 0.;
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_boundary(int p, fields_t res, int mb, int me)
{
  // lo
  for (int d = 0; d < 3; d++) {
    if (ppsc->patch[p].off[d] == 0) {
      if (ppsc->domain.bnd_part_lo[d] == BND_PART_REFLECTING) {
	add_ghosts_reflecting_lo(p, res, d, mb, me);
      }
    }
  }
  // hi
  for (int d = 0; d < 3; d++) {
    if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
      if (ppsc->domain.bnd_part_hi[d] == BND_PART_REFLECTING) {
	add_ghosts_reflecting_hi(p, res, d, mb, me);
      }
    }
  }
}

static void
run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	struct psc_mparticles *mprts_base, struct psc_mfields *mres,
	void (*do_run)(int p, fields_t flds, mparticles_t::patch_t prts))
{
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  mfields_t mf_res(mres);
  
  for (int p = 0; p < mprts.n_patches(); p++) {
    mf_res[p].zero();
    do_run(p, mf_res[p], mprts[p].range());
    add_ghosts_boundary(p, mf_res[p], 0, mres->nr_fields);
  }

  mprts.put_as(mprts_base, MP_DONT_COPY);
}

// ======================================================================
// n

static void
do_n_run(int p, fields_t flds, mparticles_t::patch_t prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int m = prt->kind();
    DEPOSIT_TO_GRID_2ND_NC(prt, flds, m, 1.f);
  }
}

static void
n_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
	  struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds, mprts_base, mres, do_n_run);
}

// ======================================================================
// rho

static void
do_rho_run(int p, fields_t flds, mparticles_t::patch_t prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int m = prt->kind();
    DEPOSIT_TO_GRID_2ND_NC(prt, flds, 0, ppsc->kinds[m].q);
  }
}

static void
rho_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
	    struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds, mprts_base, mres, do_rho_run);
}

