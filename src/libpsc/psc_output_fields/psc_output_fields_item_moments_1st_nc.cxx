
#include "psc_output_fields_item_private.h"

#include <math.h>

using real_t = mparticles_t::real_t;

#include "common_moments.cxx"

static void
run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	struct psc_mparticles *mprts_base, struct psc_mfields *mres_base,
	void (*do_run)(int p, fields_t flds, mparticles_t::patch_t& prts))
{
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  mfields_t mf_res = mres_base->get_as<mfields_t>(0, 0);
  for (int p = 0; p < mprts.n_patches(); p++) {
    mf_res[p].zero();
    do_run(p, mf_res[p], mprts[p]);
  }

  mf_res.put_as(mres_base, 0, mres_base->nr_fields);
  mprts.put_as(mprts_base, MP_DONT_COPY);
}

// ======================================================================
// n

static void
do_n_run(int p, fields_t flds, mparticles_t::patch_t& prts)
{
  const Grid_t& grid = ppsc->grid;
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int m = prt->kind();
    DEPOSIT_TO_GRID_1ST_NC(prt, flds, m, 1.f);
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
do_rho_run(int p, fields_t flds, mparticles_t::patch_t& prts)
{
  const Grid_t& grid = ppsc->grid;
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int m = prt->kind();
    DEPOSIT_TO_GRID_1ST_NC(prt, flds, 0, ppsc->kinds[m].q);
  }
}

static void
rho_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
	  struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds, mprts_base, mres, do_rho_run);
}

// ======================================================================
// v

static void
do_v_run(int p, fields_t flds, mparticles_t::patch_t& prts)
{
  const Grid_t& grid = ppsc->grid;
  struct psc_patch *patch = &ppsc->patch[p];
  real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
    particle_t *prt = &*prt_iter;
    int mm = prt->kind() * 3;

    real_t vxi[3];
    particle_calc_vxi(prt, vxi);

    for (int m = 0; m < 3; m++) {
      DEPOSIT_TO_GRID_1ST_NC(prt, flds, mm + m, vxi[m]);
    }
  }
}

static void
v_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds,
	  struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  run_all(item, mflds, mprts_base, mres, do_v_run);
}



#define MAKE_OP1(TYPE, NAME, FNAME, RUN)				\
  struct psc_output_fields_item_ops_##NAME##TYPE : psc_output_fields_item_ops { \
  psc_output_fields_item_ops_##NAME##TYPE() {					\
    name               = #NAME #TYPE;					\
    nr_comp	       = 1;						\
    fld_names[0]       = FNAME;						\
    run_all            = RUN;						\
    flags              = POFI_ADD_GHOSTS | POFI_BY_KIND;		\
  }									\
  } psc_output_fields_item_##NAME##TYPE##_ops;

#define MAKE_OP1a(TYPE, NAME, FNAME, RUN)				\
  struct psc_output_fields_item_ops_##NAME##TYPE : psc_output_fields_item_ops { \
  psc_output_fields_item_ops_##NAME##TYPE() {					\
    name               = #NAME #TYPE;					\
    nr_comp	       = 3;						\
    fld_names[0]       = FNAME;						\
    run_all            = RUN;						\
    flags              = POFI_ADD_GHOSTS;				\
  }									\
  } psc_output_fields_item_##NAME##TYPE##_ops;

#define MAKE_OP3(TYPE, NAME, FNAMEX, FNAMEY, FNAMEZ, RUN)		\
  struct psc_output_fields_item_ops_##NAME##TYPE : psc_output_fields_item_ops { \
  psc_output_fields_item_ops_##NAME##TYPE() {					\
    name               = #NAME #TYPE;					\
    nr_comp	       = 6;						\
    fld_names[0]       = FNAMEX;					\
    fld_names[1]       = FNAMEY;					\
    fld_names[2]       = FNAMEZ;					\
    run_all            = RUN;						\
    flags              = POFI_ADD_GHOSTS | POFI_BY_KIND;		\
  }									\
  } psc_output_fields_item_##NAME##TYPE##_ops;

#define MAKE_POFI_OPS(TYPE)						\
  MAKE_OP1(TYPE, n_1st_nc_, "n_nc", n_run_all)				\
  MAKE_OP1a(TYPE, rho_1st_nc_, "rho_nc", rho_run_all)			\
  MAKE_OP3(TYPE, v_1st_nc_, "vx_nc", "vy_nc", "vz_nc", v_run_all)	\

