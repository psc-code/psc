
#pragma once

#include "fields.hxx"
#include "fields_item.hxx"
#include "checks.hxx"
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"

#include <mrc_io.h>
#include <gtensor/reductions.h>

#ifdef PSC_HAVE_ADIOS2

#include "writer_adios2.hxx"
using WriterDefault = WriterADIOS2;

#else

#include "writer_mrc.hxx"
using WriterDefault = WriterMRC;

#endif

namespace psc
{

namespace helper
{

template <typename E1, typename E2>
void print_diff(const E1& e1, const E2& e2, double eps)
{
  auto&& h_e1 = gt::host_mirror(e1);
  auto&& h_e2 = gt::host_mirror(e2);
  gt::copy(gt::eval(e1), h_e1);
  gt::copy(gt::eval(e2), h_e2);
  gt::launch_host<5>(h_e1.shape(), [=](int i, int j, int k, int m, int p) {
    auto val_e1 = h_e1(i, j, k, m, p);
    auto val_e2 = h_e2(i, j, k, m, p);
    if (std::abs(val_e1 + val_e2) > eps) {
      mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", p, i, j, k, val_e1, val_e2,
              val_e1 - val_e2);
    }
  });
}

template <typename E1, typename E2>
void print_diff_3d(const E1& e1, const E2& e2, double eps)
{
  auto&& h_e1 = gt::host_mirror(e1);
  auto&& h_e2 = gt::host_mirror(e2);
  gt::copy(gt::eval(e1), h_e1);
  gt::copy(gt::eval(e2), h_e2);
  gt::launch_host<3>(h_e1.shape(), [=](int i, int j, int k) {
    auto val_e1 = h_e1(i, j, k);
    auto val_e2 = h_e2(i, j, k);
    if (std::abs(val_e1 + val_e2) > eps) {
      mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", i, j, k, val_e1, val_e2,
              val_e1 - val_e2);
    }
  });
}

} // namespace helper

namespace checks
{

// ======================================================================
// psc::checks::continuity: Charge Continuity

template <typename S, typename ITEM_RHO>
class continuity : ChecksParams
{
public:
  using storage_type = S;
  using Item_rho = ITEM_RHO;

  continuity(const ChecksParams& params) : ChecksParams(params) {}

  // ----------------------------------------------------------------------
  // before_particle_push

  template <typename Mparticles>
  void before_particle_push(const Mparticles& mprts)
  {
    const auto& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    auto item_rho = Item_rho{grid};
    rho_m_ = psc::mflds::interior(grid, item_rho(mprts));
  }

  // ----------------------------------------------------------------------
  // after_particle_push

  template <typename Mparticles, typename MfieldsState>
  void after_particle_push(const Mparticles& mprts, MfieldsState& mflds)
  {
    const Grid_t& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    auto item_rho = Item_rho{grid};
    auto item_divj = Item_divj<MfieldsState>{};

    auto rho_p = psc::mflds::interior(grid, item_rho(mprts));
    auto divj = psc::mflds::interior(grid, item_divj(mflds));
    auto d_rho = rho_p - rho_m_;
    auto dt_divj = grid.dt * divj;

    double local_err = gt::norm_linf(d_rho + dt_divj);
    // find global max
    double max_err;
    MPI_Allreduce(&local_err, &max_err, 1, MPI_DOUBLE, MPI_MAX, grid.comm());

    if (max_err >= continuity_threshold) {
      psc::helper::print_diff(d_rho, -dt_divj, continuity_threshold);
    }

    if (continuity_verbose || max_err >= continuity_threshold) {
      mpi_printf(grid.comm(), "continuity: max_err = %g (thres %g)\n", max_err,
                 continuity_threshold);
    }

    if (continuity_dump_always || max_err >= continuity_threshold) {
      if (!writer_) {
        writer_.open("continuity");
      }
      writer_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_.write(dt_divj, grid, "dt_divj", {"dt_divj"});
      writer_.write(d_rho, grid, "d_rho", {"d_rho"});
      writer_.end_step();
      MPI_Barrier(grid.comm());
    }

    assert(max_err < continuity_threshold);
  }

private:
  storage_type rho_m_;
  WriterDefault writer_;
};

// ======================================================================
// psc::checks::gauss: Gauss's Law div E = rho

template <typename S, typename ITEM_RHO>
class gauss : ChecksParams
{
public:
  using storage_type = S;
  using Item_rho = ITEM_RHO;

  gauss(const ChecksParams& params) : ChecksParams(params) {}

  // ----------------------------------------------------------------------
  // operator()

  template <typename Mparticles, typename MfieldsState>
  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    const auto& grid = mprts.grid();
    if (gauss_every_step <= 0 || grid.timestep() % gauss_every_step != 0) {
      return;
    }

    auto item_rho = Item_rho{grid};
    auto item_dive = Item_dive<MfieldsState>{};
    auto rho = psc::mflds::interior(grid, item_rho(mprts));
    auto dive = psc::mflds::interior(grid, item_dive(mflds));

    double max_err = 0.;
    for (int p = 0; p < grid.n_patches(); p++) {
      int l[3] = {0, 0, 0}, r[3] = {0, 0, 0};
      for (int d = 0; d < 3; d++) {
        if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
            grid.atBoundaryLo(p, d)) {
          l[d] = 1;
        }
      }

      auto patch_rho =
        rho.view(_s(l[0], -r[0]), _s(l[1], -r[1]), _s(l[2], -r[2]), 0, p);
      auto patch_dive =
        dive.view(_s(l[0], -r[0]), _s(l[1], -r[1]), _s(l[2], -r[2]), 0, p);

      double patch_err = gt::norm_linf(patch_dive - patch_rho);
      max_err = std::max(max_err, patch_err);

      if (patch_err > gauss_threshold) {
        psc::helper::print_diff_3d(patch_rho, patch_dive, gauss_threshold);
      }
    }

    // find global max
    double tmp = max_err;
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, grid.comm());

    if (gauss_verbose || max_err >= gauss_threshold) {
      mpi_printf(grid.comm(), "gauss: max_err = %g (thres %g)\n", max_err,
                 gauss_threshold);
    }

    if (gauss_dump_always || max_err >= gauss_threshold) {
      if (!writer_) {
        writer_.open("gauss");
      }
      writer_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_.write(rho, grid, "rho", {"rho"});
      writer_.write(dive, grid, "dive", {"dive"});
      writer_.end_step();
    }

    assert(max_err < gauss_threshold);
  }

private:
  WriterDefault writer_;
};

} // namespace checks
} // namespace psc

struct checks_order_1st
{
  template <typename S, typename D>
  using Moment_rho_nc = Moment_rho_1st_nc<S, D>;
};

struct checks_order_2nd
{
  template <typename S, typename D>
  using Moment_rho_nc = Moment_rho_2nd_nc<S, D>;
};

template <typename MP, typename S, typename ITEM_RHO>
class ChecksCommon : public ChecksParams
{
public:
  using Mparticles = MP;
  using storage_type = S;
  using item_rho_type = ITEM_RHO;

  ChecksCommon(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams{params}, continuity_{params}, gauss_{params}
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
  psc::checks::continuity<storage_type, item_rho_type> continuity_;
  psc::checks::gauss<storage_type, item_rho_type> gauss_;
};

template <typename MP, typename MF, typename ORDER, typename D>
using Checks_ =
  ChecksCommon<MP, typename MF::Storage,
               typename ORDER::template Moment_rho_nc<typename MF::Storage, D>>;
