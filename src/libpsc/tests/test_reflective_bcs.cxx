#include "gtest/gtest.h"

#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "add_ghosts_reflecting.hxx"
#include "../psc_config.hxx"

// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

// ----------------------------------------------------------------------

using BgkMfields = PscConfig::Mfields;
using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;

// ======================================================================
// Global parameters

namespace
{
PscParams psc_params;
} // namespace

// ======================================================================
// setupGrid

Grid_t* setupGrid()
{
  auto domain = Grid_t::Domain{{1, 8, 2},          // n grid points
                               {10.0, 80.0, 20.0}, // physical lengths
                               {0, 0, 0},          // location of lower corner
                               {1, 1, 1}};         // n patches

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 1.0, "i"};

  // --- generic setup
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 1;

  double dt = psc_params.cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  Int3 ibn = {2, 2, 2};
  if (Dim::InvarX::value) {
    ibn[0] = 0;
  }
  if (Dim::InvarY::value) {
    ibn[1] = 0;
  }
  if (Dim::InvarZ::value) {
    ibn[2] = 0;
  }

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

// ======================================================================
// Integration test (pun not intended)

TEST(ReflectiveBcsTest, Integration)
{
  // ----------------------------------------------------------------------
  // setup

  psc_params.nmax = 100;
  psc_params.stats_every = 1;
  psc_params.cfl = .75;

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  ChecksParams checks_params{};
  checks_params.continuity.check_interval = 1;
  checks_params.gauss.check_interval = 1;
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  Balance balance{.1};
  Collision collision{grid, 0, 0.1};
  Marder marder(grid, 0.9, 3, false);

  OutputFields<MfieldsState, Mparticles, Dim> outf{grid, {}};
  OutputParticles outp{grid, {}};
  DiagEnergies oute{grid.comm(), 0};
  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  auto psc =
    makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts, balance,
                                 collision, checks, marder, diagnostics);

  // ----------------------------------------------------------------------
  // set up initial conditions

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  {
    auto injector = mprts.injector();
    auto inj = injector[p];
    double vy = .1;
    double y_center = grid.domain.length[1] / 2.0;
    // inject 2 particles at same location to satisfy Gauss' law at t=0
    inj({{0, y_center, 5}, {0, -vy, 0}, 0, 0}); // electron
    inj({{0, y_center, 5}, {0, vy, 0}, 1, 0});  // positron
  }

  // ----------------------------------------------------------------------
  // run the simulation

  auto prts = mprts.accessor()[p];

  EXPECT_EQ(prts.size(), 2);
  EXPECT_LT(prts[0].u()[1], 0.0);
  EXPECT_GT(prts[1].u()[1], 0.0);

  for (; grid.timestep_ < psc_params.nmax; grid.timestep_++) {
    psc.step();
    ASSERT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    ASSERT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  EXPECT_EQ(prts.size(), 2);
  EXPECT_GT(prts[0].u()[1], 0.0);
  EXPECT_LT(prts[1].u()[1], 0.0);
}

// ======================================================================
// Unit tests

int ravel_idx(Int3 idx, Int3 shape)
{
  return idx[0] + shape[0] * (idx[1] + shape[1] * (idx[2]));
}

template <typename FE>
void init_mres(FE& mres_gt, Int3 shape, int p)
{
  for (int i = 0; i < shape[0]; i++) {
    for (int j = 0; j < shape[1]; j++) {
      for (int k = 0; k < shape[2]; k++) {
        int cell_id = ravel_idx({i, j, k}, shape);
        mres_gt(i, j, k, 0, p) = cell_id;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingHighY)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 1;
  add_ghosts_reflecting_hi(ldims, mres, ib, p, dim, 0, 1);

  int bx = -ib[0];
  int by = -ib[1];
  int bz = -ib[2];
  for (int x = -bx; x < ldims[0] + bx; x++) {
    for (int y = ldims[1] - by; y < ldims[1]; y++) {
      for (int z = -bz; z < ldims[2] + bz; z++) {
        int cell_id = ravel_idx({bx + x, by + y, bz + z}, shape);

        int reflected_y = 2 * ldims[1] - y - 1;
        int reflected_cell_id =
          ravel_idx({bx + x, by + reflected_y, bz + z}, shape);
        EXPECT_EQ(mres(bx + x, by + y, bz + z, 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingLowY)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 1;
  add_ghosts_reflecting_lo(ldims, mres, ib, p, dim, 0, 1);

  int bx = -ib[0];
  int by = -ib[1];
  int bz = -ib[2];
  for (int x = -bx; x < ldims[0] + bx; x++) {
    for (int y = 0; y < by; y++) {
      for (int z = -bz; z < ldims[2] + bz; z++) {
        int cell_id = ravel_idx({bx + x, by + y, bz + z}, shape);

        int reflected_y = -y - 1;
        int reflected_cell_id =
          ravel_idx({bx + x, by + reflected_y, bz + z}, shape);
        EXPECT_EQ(mres(bx + x, by + y, bz + z, 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingHighZ)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 2;
  add_ghosts_reflecting_hi(ldims, mres, ib, p, dim, 0, 1);

  int bx = -ib[0];
  int by = -ib[1];
  int bz = -ib[2];
  for (int x = -bx; x < ldims[0] + bx; x++) {
    for (int y = -by; y < ldims[1] + by; y++) {
      for (int z = ldims[2] - bz; z < ldims[2]; z++) {
        int cell_id = ravel_idx({bx + x, by + y, bz + z}, shape);

        int reflected_z = 2 * ldims[2] - z - 1;
        int reflected_cell_id =
          ravel_idx({bx + x, by + y, bz + reflected_z}, shape);
        EXPECT_EQ(mres(bx + x, by + y, bz + z, 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingLowZ)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 2;
  add_ghosts_reflecting_lo(ldims, mres, ib, p, dim, 0, 1);

  int bx = -ib[0];
  int by = -ib[1];
  int bz = -ib[2];
  for (int x = -bx; x < ldims[0] + bx; x++) {
    for (int y = -by; y < ldims[1] + by; y++) {
      for (int z = 0; z < bz; z++) {
        int cell_id = ravel_idx({bx + x, by + y, bz + z}, shape);

        int reflected_z = -z - 1;
        int reflected_cell_id =
          ravel_idx({bx + x, by + y, bz + reflected_z}, shape);
        EXPECT_EQ(mres(bx + x, by + y, bz + z, 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
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
