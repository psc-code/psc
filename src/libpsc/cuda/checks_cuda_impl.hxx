
#pragma once

#include "checks.hxx"

#include "fields_item_moments_1st_cuda.hxx"

template <typename MP, typename D>
struct ChecksCuda : ChecksParams
{
  using Mparticles = MP;
  using dim_t = D;
  using storage_type = MfieldsCuda::Storage;
  using Moment_t = Moment_rho_1st_nc_cuda<D>;

  ChecksCuda(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params), continuity_{params}, gauss_{params}
  {}

  void continuity_before_particle_push(Mparticles& mprts)
  {
    continuity_.before_particle_push(mprts);
  }

  template <typename MfieldsState>
  void continuity_after_particle_push(Mparticles& mprts, MfieldsState& mflds)
  {
    continuity_.after_particle_push(mprts, mflds);
  }

  template <typename MfieldsState>
  void gauss(Mparticles& mprts, MfieldsState& mflds)
  {
    gauss_(mprts, mflds);
  }

private:
  psc::checks::continuity<storage_type, Moment_t> continuity_;
  psc::checks::gauss<storage_type, Moment_t> gauss_;
};
