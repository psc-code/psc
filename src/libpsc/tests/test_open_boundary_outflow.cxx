#include <gtest/gtest.h>

#include "test_common.hxx"
#include "open_boundary_test_grid.hxx"

#include "psc.hxx"
#include "../psc_config.hxx"

// Particles leaving through open boundaries, with no injectors: their virtual
// charge must end up beyond the boundary nodes, so Gauss's law holds on every
// node of the domain.

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;

// 8 x 8 cells, open in y and z, unless `open` is false (then periodic
// everywhere)
Grid_t* setupGrid(bool open = true)
{
  return makeOpenBoundaryTestGrid({1, 8, 8}, {false, open, open}, 1.);
}

ChecksParams makeChecksParams()
{
  ChecksParams checks_params{};
  checks_params.gauss.check_interval = 1;
  return checks_params;
}

// Places a neutral pair at `x`, gives each particle its momentum, steps once,
// and checks Gauss's law.
template <typename Check>
void run_pair(Vec3<double> x, Vec3<double> u_e, Vec3<double> u_i, Check check)
{
  PscParams psc_params;
  psc_params.nmax = 1;

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};
  Checks checks{grid, MPI_COMM_WORLD, makeChecksParams()};
  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);
  {
    auto injector = mprts.injector();
    auto inj = injector[0];
    inj({x, u_e, 1, KIND_ELECTRON});
    inj({x, u_i, 1, KIND_ION});
  }

  psc.pre_first_step();
  psc.step();
  EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  check(grid, mprts, mflds);
}

TEST(OpenBoundaryOutflowTest, Lower)
{
  // the electron stops at y = -.594 and is deleted
  run_pair({.5, .3, 4.5}, {0, -2, 0}, {0, 0, 0},
           [](const Grid_t& grid, Mparticles& mprts, MfieldsState& mflds) {
             ASSERT_EQ(mprts.size(), 1);
           });
}

TEST(OpenBoundaryOutflowTest, Upper)
{
  // both leave, to different depths in the ghost cell. No real charge is
  // left, so no virtual charge may remain on the boundary node either: E_y
  // vanishes on both sides of it (the upper boundary nodes aren't interior)
  run_pair({.5, 7.7, 4.5}, {0, 2, 0}, {0, 1, 0},
           [](const Grid_t& grid, Mparticles& mprts, MfieldsState& mflds) {
             ASSERT_EQ(mprts.size(), 0);
             auto F = make_Fields3d<Dim>(mflds[0]);
             int ny = grid.ldims[1];
             for (int k = 0; k < grid.ldims[2]; k++) {
               EXPECT_NEAR(F(EY, 0, ny - 1, k), 0., 1e-10);
               EXPECT_NEAR(F(EY, 0, ny, k), 0., 1e-10);
             }
           });
}

TEST(OpenBoundaryOutflowTest, Corner)
{
  // the electron stops at (y, z) = (-.367, -.367), beyond two boundaries
  run_pair({.5, .3, .3}, {0, -2, -2}, {0, 0, 0},
           [](const Grid_t& grid, Mparticles& mprts, MfieldsState& mflds) {
             ASSERT_EQ(mprts.size(), 1);
           });
}

TEST(OpenBoundaryOutflowTest, MarderAtOpenBoundaries)
{
  // Marder corrects open boundary nodes; with nothing to correct, it must leave
  // the fields (ghosts included) finite and Gauss's law intact
  PscParams psc_params;
  psc_params.nmax = 3;
  psc_params.marder_interval = 1;

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};
  Checks checks{grid, MPI_COMM_WORLD, makeChecksParams()};
  Balance balance{.1};
  Collision collision{grid, 0, 0.1};
  PscConfig::Marder marder(grid, .9, 3, false);

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);
  psc.add_gauss_corrector(&marder);
  {
    auto injector = mprts.injector();
    auto inj = injector[0];
    inj({{.5, .3, 4.5}, {0, -2, 0}, 1, KIND_ELECTRON});
    inj({{.5, .3, 4.5}, {0, 0, 0}, 1, KIND_ION});
  }

  psc.pre_first_step();
  for (; grid.timestep_ < psc_params.nmax;) {
    psc.step();
    EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  auto E = mflds.storage().view(_all, _all, _all, _s(EX, EX + 3), _all);
  EXPECT_LT(gt::norm_linf(E), 1.);
}

TEST(OpenBoundaryOutflowTest, SecondOrderPscWithPeriodicBoundaries)
{
  // the 2nd-order pusher has no open boundary outflow, which mustn't keep it
  // from building or running without open boundaries
  using PscConfig2nd = PscConfig2ndDouble<Dim>;

  PscParams psc_params;
  psc_params.nmax = 1;

  auto grid_ptr = setupGrid(false);
  auto& grid = *grid_ptr;
  PscConfig2nd::MfieldsState mflds{grid};
  PscConfig2nd::Mparticles mprts{grid};
  PscConfig2nd::Checks checks{grid, MPI_COMM_WORLD, ChecksParams{}};
  PscConfig2nd::Balance balance{.1};
  PscConfig2nd::Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig2nd>(psc_params, grid, mflds, mprts,
                                             balance, collision, checks);
  {
    auto injector = mprts.injector();
    auto inj = injector[0];
    inj({{.5, 4.5, 4.5}, {0, 1, 0}, 1, KIND_ELECTRON});
    inj({{.5, 4.5, 4.5}, {0, 0, 0}, 1, KIND_ION});
  }

  psc.pre_first_step();
  psc.step();
  EXPECT_EQ(mprts.size(), 2);
}

TEST(OpenBoundaryOutflowTest, SinglePrecision)
{
  // an electron in the lower y ghost cell deposits its exit current
  using PscConfigSingle = PscConfig1vbecSingle<Dim>;
  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  PscConfigSingle::MfieldsState mflds{grid};
  PscConfigSingle::Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0].inject_local_unchecked(
      {{.5, -.3, 4.5}, {0, -2, 0}, 1, KIND_ELECTRON});
  }

  PscConfigSingle::PushParticles::OpenBoundaryOutflow_t outflow;
  outflow.inject(mprts, mflds);

  auto F = make_Fields3d<Dim>(mflds[0]);
  EXPECT_NE(F(JXI + 1, 0, -1, 4), 0.f);
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
