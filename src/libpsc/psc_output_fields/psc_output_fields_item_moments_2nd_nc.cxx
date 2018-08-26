
#include "psc_output_fields_item_private.h"
#include "fields.hxx"

#include <math.h>

#include "common_moments.cxx"

// ======================================================================
// n

template<typename MP, typename MF>
struct Moment_n_2nd_nc
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::patch_t;
  using fields_t = typename Mfields::fields_t;  
  
  constexpr static char const* name = "n_2nd_nc_double";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "n" }; }
  constexpr static int flags = POFI_BY_KIND;
  
  static void run(fields_t flds, particles_t& prts)
  {
    const auto& grid = prts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      auto *prt = &*prt_iter;
      int m = prt->kind;
      DEPOSIT_TO_GRID_2ND_NC(prt, flds, m, 1.f);
    }
  }
};

// ======================================================================
// rho

template<typename MP, typename MF>
struct Moment_rho_2nd_nc
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::patch_t;
  using fields_t = typename Mfields::fields_t;  
  
  constexpr static char const* name = "rho_2nd_nc";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "rho" }; }
  constexpr static int flags = 0;
  
  static void run(fields_t flds, particles_t& prts)
  {
    const auto& grid = flds.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      auto *prt = &*prt_iter;
      DEPOSIT_TO_GRID_2ND_NC(prt, flds, 0, prts.prt_qni(*prt));
    }
  }
};

