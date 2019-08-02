
#pragma once

#include <math.h>

#include "common_moments.cxx"

// ======================================================================
// n

template<typename MP, typename MF>
struct Moment_n_1st_nc
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  
  constexpr static char const* name = "n_1st_nc";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return { "n" }; }
  constexpr static int flags = POFI_BY_KIND;
  
  static void run(Mfields& mflds, Mparticles& mprts)
  {
    const auto& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];

    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      for (auto prt: accessor[p]) {
	int m = prt.kind();
	DEPOSIT_TO_GRID_1ST_NC(prt, flds, m, 1.f);
      }
    }
  }
};

// ======================================================================
// rho

template<typename MP, typename MF>
struct Moment_rho_1st_nc : ItemMomentCRTP<Moment_rho_1st_nc<MP, MF>, MF>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc<MP, MF>, MF>;
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  
  constexpr static char const* name = "rho_1st_nc";
  static int n_comps(const Grid_t& grid) { return 1; }
  static std::vector<std::string> fld_names() { return { "rho" }; }
  constexpr static int flags = 0;

  explicit Moment_rho_1st_nc(const Mparticles& mprts) : Base{mprts.grid()}
  {
    const auto& grid = mprts.grid();
    real_t fnqs = grid.norm.fnqs;
    real_t dxi = 1.f / grid.domain.dx[0], dyi = 1.f / grid.domain.dx[1], dzi = 1.f / grid.domain.dx[2];
    
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto res = Base::mres_[p];
      for (auto prt: accessor[p]) {
	DEPOSIT_TO_GRID_1ST_NC(prt, res, 0, prt.q());
      }
    }
    Base::bnd_.add_ghosts(Base::mres_);
  }
};

// ======================================================================
// v

template<typename MP, typename MF>
struct Moment_v_1st_nc
{
  using Mparticles = MP;
  using Mfields = MF;
  using real_t = typename Mparticles::real_t;
  using particles_t = typename Mparticles::Patch;
  
  constexpr static char const* name = "v_1st_nc";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return { "vx", "vy", "vz" }; }
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
	  DEPOSIT_TO_GRID_1ST_NC(prt, flds, mm + m, vxi[m]);
	}
      }
    }
  }
};

#define MAKE_POFI_OPS(MP, MF, TYPE)					\
  FieldsItemMomentOps<Moment_n_1st_nc<MP, MF>> psc_output_fields_item_n_1st_nc_##TYPE##_ops; \
  FieldsItemMomentOps<Moment_rho_1st_nc<MP, MF>> psc_output_fields_item_rho_1st_nc_##TYPE##_ops; \
  FieldsItemMomentOps<Moment_v_1st_nc<MP, MF>> psc_output_fields_item_v_1st_nc_##TYPE##_ops; \

