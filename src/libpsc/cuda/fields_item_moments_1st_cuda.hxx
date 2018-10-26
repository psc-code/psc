
#pragma once

#include "psc_particles_cuda.h"
#include "fields_item.hxx"
#include "bnd_cuda_impl.hxx"
#include "cuda_moments.cuh"

// ======================================================================
// Moment_rho_1st_nc_cuda

template<typename _Mparticles, typename dim>
struct Moment_rho_1st_nc_cuda : ItemMomentCRTP<Moment_rho_1st_nc_cuda<_Mparticles, dim>, MfieldsCuda>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc_cuda, MfieldsCuda>;
  using Mparticles = _Mparticles;
  using Mfields = MfieldsCuda;
  using Bnd = BndCuda<Mfields>;
  
  constexpr static const char* name = "rho_1st_nc";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "rho_nc_cuda" }; } // FIXME
  constexpr static int flags = 0;

  Moment_rho_1st_nc_cuda(const Grid_t& grid, MPI_Comm comm)
    : Base(grid, comm),
      bnd_{grid, grid.ibn}
  {}

  void run(Mparticles& mprts)
  {
    Mfields& mres  = this->mres_;
    auto& cmprts = *mprts.cmprts();
    cuda_mfields *cmres = mres.cmflds();
    
    mres.zero();
    CudaMoments1stNcRho<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, cmres);
    bnd_.add_ghosts(mres, 0, mres.n_comps());
  }

private:
  Bnd bnd_;
};

// ======================================================================
// n_1st_cuda

template<typename _Mparticles, typename dim>
struct Moment_n_1st_cuda : ItemMomentCRTP<Moment_n_1st_cuda<_Mparticles, dim>, MfieldsCuda>
{
  using Base = ItemMomentCRTP<Moment_n_1st_cuda, MfieldsCuda>;
  using Mparticles = _Mparticles;
  using Mfields = MfieldsCuda;
  using Bnd = BndCuda<Mfields>;
  
  constexpr static const char* name = "n_1st";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "n_1st_cuda" }; }
  constexpr static int flags = POFI_BY_KIND;

  Moment_n_1st_cuda(const Grid_t& grid, MPI_Comm comm)
    : Base(grid, comm),
      bnd_{grid, grid.ibn}
  {}

  void run(Mparticles& mprts)
  {
    Mfields& mres = this->mres_;
    auto& cmprts = *mprts.cmprts();
    cuda_mfields *cmres = mres.cmflds();
    
    mres.zero();
    CudaMoments1stNcN<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, cmres);
    bnd_.add_ghosts(mres, 0, mres.n_comps());
  }

private:
  Bnd bnd_;
};

