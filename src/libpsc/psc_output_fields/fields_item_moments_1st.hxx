
#pragma once

#include <cmath>

#include "common_moments.cxx"

// ======================================================================
// n_1st

template<typename MP, typename MF>
struct Moment_n_1st
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  using fields_t = typename Mfields::fields_t;  
  
  constexpr static char const* name = "n_1st";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "n" }; }
  constexpr static int flags = POFI_BY_KIND;
  
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const Grid_t& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];

    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt: accessor[p]) {
	int m = prt.kind();
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, m, 1.f);
      }
    }
  }
};

// ======================================================================
// v_1st

template<typename MP, typename MF>
struct Moment_v_1st
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  using fields_t = typename Mfields::fields_t;  
  
  constexpr static char const* name = "v_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "vx", "vy", "vz" }; }
  constexpr static int flags = POFI_BY_KIND;

  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const Grid_t& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];
    
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt: accessor[p]) {
	int mm = prt.kind() * 3;
      
	real_t vxi[3];
	particle_calc_vxi(prt, vxi);
	
	for (int m = 0; m < 3; m++) {
	  DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, vxi[m]);
	}
      }
    }
  }
};

// ======================================================================
// p_1st

template<typename MP, typename MF>
struct Moment_p_1st
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  using fields_t = typename Mfields::fields_t;  
  
  constexpr static char const* name = "p_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "px", "py", "pz" }; }
  constexpr static int flags = POFI_BY_KIND;

  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const Grid_t& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];
    
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt: accessor[p]) {
	int mm = prt.kind() * 3;
	auto pxi = prt.u();
	
	for (int m = 0; m < 3; m++) {
	  DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, prt.m() * pxi[m]);
	}
      }
    }
  }
};

// ======================================================================
// vv_1st

template<typename MP, typename MF>
struct Moment_vv_1st
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  using fields_t = typename Mfields::fields_t;  
  
  constexpr static char const* name = "vv_1st";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "vxvx", "vyvy", "vzvz" }; }
  constexpr static int flags = POFI_BY_KIND;

  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const Grid_t& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];
    
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt: accessor[p]) {
	int mm = prt.kind() * 3;
	
	real_t vxi[3];
	particle_calc_vxi(prt, vxi);
	
	for (int m = 0; m < 3; m++) {
	  DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + m, vxi[m] * vxi[m]);
	}
      }
    }
  }
};

// ======================================================================
// T_1st

template<typename MP, typename MF>
struct Moment_T_1st
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  using fields_t = typename Mfields::fields_t;  
  
  constexpr static char const* name = "T_1st";
  constexpr static int n_comps = 6;
  constexpr static fld_names_t fld_names() { return { "Txx", "Tyy", "Tzz", "Txy", "Txz", "Tyz" }; }
  constexpr static int flags = POFI_BY_KIND;

  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const Grid_t& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];
    
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt: accessor[p]) {
	int mm = prt.kind() * 6;
	
	real_t vxi[3];
	particle_calc_vxi(prt, vxi);
	auto pxi = prt.u();
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 0, prt.m() * pxi[0] * vxi[0]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 1, prt.m() * pxi[1] * vxi[1]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 2, prt.m() * pxi[2] * vxi[2]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 3, prt.m() * pxi[0] * vxi[1]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4, prt.m() * pxi[0] * vxi[2]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 5, prt.m() * pxi[1] * vxi[2]);
      }
    }
  }
};

// ======================================================================
// Tvv_1st

template<typename MP, typename MF>
struct Moment_Tvv_1st
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  using fields_t = typename Mfields::fields_t;  
  
  constexpr static char const* name = "Tvv_1st";
  constexpr static int n_comps = 6;
  constexpr static fld_names_t fld_names() { return { "vxvx", "vyvy", "vzvz", "vxvy", "vxvz", "vyvz" }; }
  constexpr static int flags = POFI_BY_KIND;

  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const Grid_t& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];
    
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt: accessor[p]) {
	int mm = prt.kind() * 6;
      
	real_t vxi[3];
	particle_calc_vxi(prt, vxi);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 0, prt.m() * vxi[0] * vxi[0]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 1, prt.m() * vxi[1] * vxi[1]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 2, prt.m() * vxi[2] * vxi[2]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 3, prt.m() * vxi[0] * vxi[1]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 4, prt.m() * vxi[0] * vxi[2]);
	DEPOSIT_TO_GRID_1ST_CC(prt, flds, mm + 5, prt.m() * vxi[1] * vxi[2]);
      }
    }
  }
};

