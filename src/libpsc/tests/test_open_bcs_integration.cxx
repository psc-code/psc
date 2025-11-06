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

  auto bc = psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_OPEN},
                          {BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_OPEN},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_OPEN},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_OPEN}};

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
// Start with 2 particles, send them flying off towards the boundary, and demand
// that:
// a) they get dropped
// b) there are no continuity or gauss errors in the fields

TEST(OpenBcsTest, IntegrationY)
{
  // ----------------------------------------------------------------------
  // setup

  psc_params.nmax = 20;
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
    double vy = 100; // fast enough to escape electric attraction
    double y_quarter = grid.domain.length[1] / 4.0;
    // inject 2 particles at same location to satisfy Gauss' law at t=0
    inj({{0, y_quarter, 5}, {0, -vy, 0}, 1, 0}); // electron
    inj({{0, y_quarter, 5}, {0, vy, 0}, 1, 1});  // positron
  }

  // ----------------------------------------------------------------------
  // run the simulation

  auto accessor = mprts.accessor();
  auto prts = accessor[p];

  ASSERT_EQ(prts.size(), 2);
  ASSERT_LT(prts[0].u()[1], 0.0);
  ASSERT_GT(prts[1].u()[1], 0.0);

  ASSERT_EQ(prts[0].m(), 1.0);
  ASSERT_EQ(prts[1].m(), 1.0);

  ASSERT_EQ(prts[0].q(), -1.0);
  ASSERT_EQ(prts[1].q(), 1.0);

  for (; grid.timestep_ < psc_params.nmax;) {
    psc.step();
    ASSERT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    ASSERT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  ASSERT_TRUE(prts.size() < 2);
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
