#include "gtest/gtest.h"

#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
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
// run

static void run(int argc, char** argv)
{
  mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

  // ----------------------------------------------------------------------
  // setup various parameters first

  psc_params.nmax = 100;
  psc_params.stats_every = 1;
  psc_params.cfl = .75;

  // ----------------------------------------------------------------------
  // Set up grid, state fields, particles

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  // ----------------------------------------------------------------------
  // Set up various objects needed to run this case

  // -- Balance
  psc_params.balance_interval = 0;
  Balance balance{.1};

  // -- Sort
  psc_params.sort_interval = 100;

  // -- Collision
  int collision_interval = 0;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.gauss.check_interval = 1;
  checks_params.continuity.check_interval = 1;

  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = -1;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output

  // -- output fields
  OutputFieldsParams outf_params{};
  outf_params.fields.pfield.out_interval = -1;
  outf_params.moments.pfield.out_interval = -1;
  OutputFields<MfieldsState, Mparticles, Dim> outf{grid, outf_params};

  // -- output particles
  OutputParticlesParams outp_params{};
  outp_params.every_step = -1;
  OutputParticles outp{grid, outp_params};

  int oute_interval = -100;
  DiagEnergies oute{grid.comm(), oute_interval};

  auto diagnostics = makeDiagnosticsDefault(outf, outp, oute);

  // ----------------------------------------------------------------------
  // set up initial conditions

  {
    auto injector = mprts.injector();
    auto inj = injector[0];
    inj({{0, 15, 5}, {0, -.001, 0}, 1, 0});
  }

  // ----------------------------------------------------------------------
  // run the simulation

  auto psc =
    makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts, balance,
                                 collision, checks, marder, diagnostics);

  for (; grid.timestep_ < psc_params.nmax; grid.timestep_++) {
    psc.step();
    ASSERT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    ASSERT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }
}

// ======================================================================
// Integration test (pun not intended)

TEST(ReflectiveBcsTest, Integration)
{
  int argc = 1;
  char** argv = nullptr;

  run(argc, argv);
}

// ======================================================================
// Unit tests

int ravel_idx(Int3 idx, Int3 shape)
{
  return idx[0] + shape[0] * (idx[1] + shape[1] * (idx[2]));
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

  for (int i = 0; i < shape[0]; i++) {
    for (int j = 0; j < shape[1]; j++) {
      for (int k = 0; k < shape[2]; k++) {
        int cell_id = ravel_idx({i, j, k}, shape);
        mres(i, j, k, 0, p) = cell_id;
      }
    }
  }

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
  // TODO
  EXPECT_TRUE(false);
}

TEST(ReflectiveBcsTest, AddGhostsReflectingHighZ)
{
  // TODO
  EXPECT_TRUE(false);
}

TEST(ReflectiveBcsTest, AddGhostsReflectingLowZ)
{
  // TODO
  EXPECT_TRUE(false);
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
