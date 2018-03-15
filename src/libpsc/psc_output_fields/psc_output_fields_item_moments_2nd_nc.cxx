
#include "psc_output_fields_item_private.h"
#include "fields.hxx"

using fields_t = mfields_t::fields_t;
using Fields = Fields3d<fields_t>;
using real_t = mparticles_t::real_t;

#include <math.h>

#include "common_moments.cxx"

using particles_t = mparticles_t::patch_t;

// ======================================================================
// n

struct Moment_n
{
  using mparticles_t = mparticles_t;
  using mfields_t = mfields_t;
  
  constexpr static char const* name = "n_2nd_nc_double";
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
      DEPOSIT_TO_GRID_2ND_NC(prt, flds, m, 1.f);
    }
  }
};

// ======================================================================
// rho

struct Moment_rho
{
  using mparticles_t = mparticles_t;
  using mfields_t = mfields_t;
  
  constexpr static char const* name = "rho_2nd_nc_double";
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
      DEPOSIT_TO_GRID_2ND_NC(prt, flds, 0, ppsc->kinds[m].q);
    }
  }
};

