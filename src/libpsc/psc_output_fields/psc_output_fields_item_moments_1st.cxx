
#include "psc_output_fields_item_private.h"
#include <psc_bnd.h>
#include <fields.hxx>
#include <bnd.hxx>
#include <fields_item.hxx>

#include <string>
#include <math.h>

#include "common_moments.cxx"

// ======================================================================
// n_1st

template<typename MP, typename MF>
struct Moment_n_1st
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using real_t = typename mparticles_t::real_t;
  using particles_t = typename mparticles_t::patch_t;
  using fields_t = typename mfields_t::fields_t;  
  
  constexpr static char const* name = "n_1st";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "n" }; }
  constexpr static int flags = POFI_ADD_GHOSTS | POFI_BY_KIND;
  
  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      auto *prt = &*prt_iter;
      int m = prt->kind();
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, m, 1.f);
    }
  }
};

// ======================================================================
// v_1st

template<typename MP, typename MF>
struct Moment_v_1st
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using real_t = typename mparticles_t::real_t;
  using particles_t = typename mparticles_t::patch_t;
  using fields_t = typename mfields_t::fields_t;  
  
  constexpr static char const* name = "v_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "vx", "vy", "vz" }; }
  constexpr static int flags = POFI_ADD_GHOSTS | POFI_BY_KIND;

  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      auto *prt = &*prt_iter;
      int mm = prt->kind() * 3;
      
      real_t vxi[3];
      particle_calc_vxi(prt, vxi);
      
      for (int m = 0; m < 3; m++) {
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, vxi[m]);
      }
    }
  }
};

// ======================================================================
// p_1st

template<typename MP, typename MF>
struct Moment_p_1st
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using real_t = typename mparticles_t::real_t;
  using particles_t = typename mparticles_t::patch_t;
  using fields_t = typename mfields_t::fields_t;  
  
  constexpr static char const* name = "p_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "px", "py", "pz" }; }
  constexpr static int flags = POFI_ADD_GHOSTS | POFI_BY_KIND;

  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      auto *prt = &*prt_iter;
      int mm = prt->kind() * 3;
      real_t *pxi = &prt->pxi;
      
      for (int m = 0; m < 3; m++) {
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, prts.prt_mni(*prt) * pxi[m]);
      }
    }
  }
};

// ======================================================================
// vv_1st

template<typename MP, typename MF>
struct Moment_vv_1st
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using real_t = typename mparticles_t::real_t;
  using particles_t = typename mparticles_t::patch_t;
  using fields_t = typename mfields_t::fields_t;  
  
  constexpr static char const* name = "vv_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "vxvx", "vyvy", "vzvz" }; }
  constexpr static int flags = POFI_ADD_GHOSTS | POFI_BY_KIND;

  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      auto *prt = &*prt_iter;
      int mm = prt->kind() * 3;
      
      real_t vxi[3];
      particle_calc_vxi(prt, vxi);
      
      for (int m = 0; m < 3; m++) {
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, vxi[m] * vxi[m]);
      }
    }
  }
};

// ======================================================================
// T_1st

template<typename MP, typename MF>
struct Moment_T_1st
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using real_t = typename mparticles_t::real_t;
  using particles_t = typename mparticles_t::patch_t;
  using fields_t = typename mfields_t::fields_t;  
  
  constexpr static char const* name = "T_1st";
  constexpr static int n_comps = 6;
  constexpr static fld_names_t fld_names() { return { "Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz" }; }
  constexpr static int flags = POFI_ADD_GHOSTS | POFI_BY_KIND;

  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      auto *prt = &*prt_iter;
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
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 0, prts.prt_mni(*prt) * px[0] * vx[0]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 1, prts.prt_mni(*prt) * px[1] * vx[1]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 2, prts.prt_mni(*prt) * px[2] * vx[2]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 3, prts.prt_mni(*prt) * px[0] * vx[1]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4, prts.prt_mni(*prt) * px[0] * vx[2]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 5, prts.prt_mni(*prt) * px[1] * vx[2]);
    }
  }
};

// ======================================================================
// Tvv_1st

template<typename MP, typename MF>
struct Moment_Tvv_1st
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using real_t = typename mparticles_t::real_t;
  using particles_t = typename mparticles_t::patch_t;
  using fields_t = typename mfields_t::fields_t;  
  
  constexpr static char const* name = "Tvv_1st";
  constexpr static int n_comps = 6;
  constexpr static fld_names_t fld_names() { return { "vxvx", "vyvy", "vzvz", "vxvy", "vxvz", "vyvz" }; }
  constexpr static int flags = POFI_ADD_GHOSTS | POFI_BY_KIND;

  static void run(fields_t flds, particles_t& prts)
  {
    const Grid_t& grid = ppsc->grid();
    real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
    real_t dxi = 1.f / grid.dx[0], dyi = 1.f / grid.dx[1], dzi = 1.f / grid.dx[2];
    
    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      auto *prt = &*prt_iter;
      int mm = prt->kind() * 6;
      
      real_t vxi[3];
      particle_calc_vxi(prt, vxi);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 0, prts.prt_mni(*prt) * vxi[0] * vxi[0]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 1, prts.prt_mni(*prt) * vxi[1] * vxi[1]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 2, prts.prt_mni(*prt) * vxi[2] * vxi[2]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 3, prts.prt_mni(*prt) * vxi[0] * vxi[1]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4, prts.prt_mni(*prt) * vxi[0] * vxi[2]);
      DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 5, prts.prt_mni(*prt) * vxi[1] * vxi[2]);
    }
  }
};

// ======================================================================

#define MAKE_OP(MP, MF, TYPE, NAME, Moment_t)				\
  FieldsItemMomentOps<Moment_t<MP, MF>> psc_output_fields_item_##NAME##TYPE##_ops;

#define MAKE_POFI_OPS(MP, MF, TYPE)					\
  MAKE_OP(MP, MF, TYPE, n_1st_  , Moment_n_1st)				\
  MAKE_OP(MP, MF, TYPE, v_1st_  , Moment_v_1st)				\
  MAKE_OP(MP, MF, TYPE, p_1st_  , Moment_p_1st)				\
  MAKE_OP(MP, MF, TYPE, vv_1st_ , Moment_vv_1st)			\
  MAKE_OP(MP, MF, TYPE, T_1st_  , Moment_T_1st)				\
  MAKE_OP(MP, MF, TYPE, Tvv_1st_, Moment_Tvv_1st)			\

