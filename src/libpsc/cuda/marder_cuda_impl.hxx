
#pragma once

#include "../libpsc/psc_push_fields/marder_impl.hxx"
#include "psc_particles_single.h"
#include "mparticles_cuda.hxx"
//#include "fields_item_dive_cuda.hxx"
#include "fields_item_moments_1st_cuda.hxx"
#include "cuda_bits.h"

#include <gtensor/reductions.h>

#include <mrc_io.h>

void cuda_marder_correct_yz(MfieldsStateCuda& mflds, MfieldsCuda& mf, int p,
                            Float3 fac, Int3 ly, Int3 ry, Int3 lz, Int3 rz);
void cuda_marder_correct_xyz(MfieldsStateCuda& mflds, MfieldsCuda& mf, int p,
                             Float3 fac, Int3 lx, Int3 rx, Int3 ly, Int3 ry,
                             Int3 lz, Int3 rz);

namespace psc
{
namespace marder
{

inline void correct(MfieldsStateCuda& mflds, MfieldsCuda& mf, float diffusion)
{
  const auto& grid = mflds.grid();

  Float3 fac;
  fac[0] = .5 * grid.dt * diffusion / grid.domain.dx[0];
  fac[1] = .5 * grid.dt * diffusion / grid.domain.dx[1];
  fac[2] = .5 * grid.dt * diffusion / grid.domain.dx[2];

  // OPT, do all patches in one kernel
  for (int p = 0; p < mflds.n_patches(); p++) {
    int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
    int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
      if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
          grid.atBoundaryLo(p, d)) {
        l_cc[d] = -1;
        l_nc[d] = -1;
      }
      if (grid.bc.fld_hi[d] == BND_FLD_CONDUCTING_WALL &&
          grid.atBoundaryHi(p, d)) {
        r_cc[d] = -1;
        r_nc[d] = 0;
      }
    }

    Int3 ldims = grid.ldims;

    Int3 lx = {l_cc[0], l_nc[1], l_nc[2]};
    Int3 rx = {r_cc[0] + ldims[0], r_nc[1] + ldims[1], r_nc[2] + ldims[2]};

    Int3 ly = {l_nc[0], l_cc[1], l_nc[2]};
    Int3 ry = {r_nc[0] + ldims[0], r_cc[1] + ldims[1], r_nc[2] + ldims[2]};

    Int3 lz = {l_nc[0], l_nc[1], l_cc[2]};
    Int3 rz = {r_nc[0] + ldims[0], r_nc[1] + ldims[1], r_cc[2] + ldims[2]};

    if (grid.isInvar(0)) {
      cuda_marder_correct_yz(mflds, mf, p, fac, ly, ry, lz, rz);
    } else {
      cuda_marder_correct_xyz(mflds, mf, p, fac, lx, rx, ly, ry, lz, rz);
    }
  }
}

} // namespace marder
} // namespace psc

template <typename BS, typename D>
struct MarderCuda
  : public MarderCommon<MparticlesCuda<BS>, MfieldsStateCuda, MfieldsCuda, D,
                        Moment_rho_1st_nc_cuda<D>, BndCuda3>
{
  using Base = MarderCommon<MparticlesCuda<BS>, MfieldsStateCuda, MfieldsCuda,
                            D, Moment_rho_1st_nc_cuda<D>, BndCuda3>;
  using Mparticles = typename Base::Mparticles;
  using MfieldsState = typename Base::MfieldsState;
  using Mfields = typename Base::Mfields;
  using dim_t = typename Base::dim_t;
  using Bnd = typename Base::Bnd;
  using real_t = typename Mfields::real_t;
  using Moment_t = typename Base::Item_rho_t;
  using Base::bnd_;
  using Base::diffusion_;
  using Base::dump_;
  using Base::io_;
  using Base::loop_;
  using Base::res_;
  using Base::rho_;

  using Base::Base;

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MfieldsStateCuda& mflds, MparticlesCuda<BS>& mprts)
  {
    const auto& grid = mprts.grid();
    static int pr;
    if (!pr) {
      pr = prof_register("marder", 1., 0, 0);
    }

    prof_start(pr);
    // need to fill ghost cells first (should be unnecessary with only variant
    // 1) FIXME
    bnd_.fill_ghosts(mflds, EX, EX + 3);

    Moment_t item_rho{grid};
    auto&& rho = psc::mflds::interior(grid, item_rho(mprts));

    for (int i = 0; i < loop_; i++) {
      Base::calc_aid_fields(mflds, rho);
      Base::print_max(res_);
      Base::correct(mflds);
      bnd_.fill_ghosts(mflds, EX, EX + 3);
    }
    prof_stop(pr);
  }
};
