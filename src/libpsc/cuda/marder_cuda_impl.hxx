
#pragma once

#include "../libpsc/psc_push_fields/marder_impl.hxx"
#include "psc_particles_single.h"
#include "mparticles_cuda.hxx"
//#include "fields_item_dive_cuda.hxx"
#include "fields_item_moments_1st_cuda.hxx"
#include "cuda_bits.h"

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

// FIXME: checkpointing won't properly restore state
// FIXME: if the subclass creates objects, it'd be cleaner to have them
// be part of the subclass

template <typename BS, typename dim>
struct MarderCuda : MarderBase
{
  using MfieldsState = MfieldsStateCuda;
  using Mfields = MfieldsCuda;
  using Mparticles = MparticlesCuda<BS>;
  using real_t = MfieldsState::real_t;
  using Moment_t = Moment_rho_1st_nc<MfieldsSingle, dim>;

  MarderCuda(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : grid_{grid},
      diffusion_{diffusion},
      loop_{loop},
      dump_{dump},
      item_rho_{grid},
      bnd_{grid, grid.ibn},
      bnd_mf_{grid, grid.ibn},
      rho_{grid, 1, grid.ibn},
      res_{grid, 1, grid.ibn}
  {
    if (dump_) {
      io_.open("marder");
    }
  }

  template <typename E>
  void calc_aid_fields(MfieldsState& mflds, const E& rho)
  {
    const auto& grid = mflds.grid();
    Item_dive<MfieldsStateCuda> item_dive(mflds);
    auto&& dive = item_dive.gt();

    if (dump_) {
      static int cnt;
      io_.begin_step(cnt, cnt); // ppsc->timestep, ppsc->timestep * ppsc->dt);
      cnt++;
      io_.write(rho, grid, "rho", {"rho"});
      io_.write(dive, grid, "dive", {"dive"});
      io_.end_step();
    }

    psc::interior(res_.gt(), res_.ib()) = dive - rho;
    bnd_mf_.fill_ghosts(res_, 0, 1);
  }

  // ----------------------------------------------------------------------
  // correct
  //
  // Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and
  // Verboncoeur, CPC, 1997)

  void correct(MfieldsState& mflds)
  {
    const Grid_t& grid = mflds._grid();
    // FIXME: how to choose diffusion parameter properly?
    float inv_sum = 0.;
    for (int d = 0; d < 3; d++) {
      if (!grid.isInvar(d)) {
        inv_sum += 1. / sqr(grid.domain.dx[d]);
      }
    }
    float diffusion_max = 1. / 2. / (.5 * grid.dt) / inv_sum;
    float diffusion = diffusion_max * diffusion_;
    psc::marder::correct(mflds, res_, diffusion);
  }

  void operator()(MfieldsStateCuda& mflds, MparticlesCuda<BS>& mprts)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("marder", 1., 0, 0);
    }

    prof_start(pr);
    // need to fill ghost cells first (should be unnecessary with only variant
    // 1) FIXME
    bnd_.fill_ghosts(mflds, EX, EX + 3);

    item_rho_(mprts);
    auto&& rho = psc::interior(item_rho_.storage(), item_rho_.ib());

    for (int i = 0; i < loop_; i++) {
      calc_aid_fields(mflds, rho);
      correct(mflds);
      bnd_.fill_ghosts(mflds, EX, EX + 3);
    }
    prof_stop(pr);
  }

  // private:
  const Grid_t& grid_;
  real_t diffusion_; //< diffusion coefficient for Marder correction
  int loop_;         //< execute this many relaxation steps in a loop
  bool dump_;        //< dump div_E, rho

  WriterMRC io_; //< for debug dumping

  BndCuda3<MfieldsState> bnd_;
  BndCuda3<Mfields> bnd_mf_;
  Mfields rho_;
  Mfields res_;

  Moment_rho_1st_nc_cuda<dim> item_rho_;
};
