
#include "psc_output_fields_item_private.h"

#include <math.h>

using real_t = mparticles_t::real_t;

#include "common_moments.cxx"

using particles_t = mparticles_t::patch_t;

// ======================================================================
// n

struct Moment_n_1st_nc
{
  using mparticles_t = mparticles_t;
  using mfields_t = mfields_t;
  
  constexpr static char const* name = "n_1st_nc";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "n" }; }
  constexpr static int flags = POFI_BY_KIND;
  
  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      particle_t *prt = &*prt_iter;
      int m = prt->kind();
      DEPOSIT_TO_GRID_1ST_NC(prt, flds, m, 1.f);
    }
  }
};

// ======================================================================
// rho

struct Moment_rho_1st_nc
{
  using mparticles_t = mparticles_t;
  using mfields_t = mfields_t;
  
  constexpr static char const* name = "rho_1st_nc";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "rho" }; }
  constexpr static int flags = 0;
  
  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      particle_t *prt = &*prt_iter;
      int m = prt->kind();
      DEPOSIT_TO_GRID_1ST_NC(prt, flds, 0, ppsc->kinds[m].q);
    }
  }
};

// ======================================================================
// v

struct Moment_v_1st_nc
{
  using mparticles_t = mparticles_t;
  using mfields_t = mfields_t;
  
  constexpr static char const* name = "v_1st_nc";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "vx", "vy", "vz" }; }
  constexpr static int flags = POFI_BY_KIND;
  
  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];

    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      particle_t *prt = &*prt_iter;
      int mm = prt->kind() * 3;
      
      real_t vxi[3];
      particle_calc_vxi(prt, vxi);
      
      for (int m = 0; m < 3; m++) {
	DEPOSIT_TO_GRID_1ST_NC(prt, flds, mm + m, vxi[m]);
      }
    }
  }
};

#define MAKE_OP(TYPE, NAME, Moment_t)					\
  FieldsItemMomentOps<Moment_t> psc_output_fields_item_##NAME##TYPE##_ops;

#define MAKE_POFI_OPS(TYPE)						\
  MAKE_OP(TYPE, n_1st_nc_  , Moment_n_1st_nc)				\
  MAKE_OP(TYPE, rho_1st_nc_, Moment_rho_1st_nc)				\
  MAKE_OP(TYPE, v_1st_nc_  , Moment_v_1st_nc)				\

