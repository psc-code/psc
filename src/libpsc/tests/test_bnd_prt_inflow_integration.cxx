#include "gtest/gtest.h"

#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "add_ghosts_reflecting.hxx"
#include "inflow.hxx"
#include "../psc_config.hxx"

// ======================================================================
// PSC configuration

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

// ----------------------------------------------------------------------

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
  norm_params.nicell = 10;

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

static double half() { return 0.5; }

// ======================================================================
// Integration test

TEST(TestSetupParticlesInflow, Integration)
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

  Double3 u = {0, 10.0, 0};
  Double3 T = {0, 0, 0};
  psc_particle_npt electron_npt = {0, 1.0, u, T};
  Inflow<Mparticles, dim_y> electron_inflow = {grid, electron_npt, *half};

  auto psc = makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts,
                                          balance, collision, checks, marder,
                                          diagnostics, electron_inflow);

  auto accessor = mprts.accessor();
  auto prts = accessor[0];

  // ----------------------------------------------------------------------
  // step 1

  psc.step();

  ASSERT_EQ(prts.size(), 20); // 20 = nicell * nz

  EXPECT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
  EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);

  // step 2

  psc.step();

  ASSERT_EQ(prts.size(), 40);

  EXPECT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
  EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  psc_init(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  psc_finalize();
  return rc;
}
