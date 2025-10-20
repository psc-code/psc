#include <gtest/gtest.h>

#include "test_common.hxx"

#include "boundary_injector.hxx"
#include "composite_injector.hxx"

#include "psc.hxx"
#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "../psc_config.hxx"

TEST(BoundaryInjectorTest, ParticleGeneratorMaxwellianTest)
{
  int kind_idx = 15;
  Grid_t::Kind kind{1.0, 1836.0, "ion"};
  ParticleGeneratorMaxwellian::Real w = 1.0;
  ParticleGeneratorMaxwellian::Real3 mean_u{0.0, 5.0, 15.0};
  ParticleGeneratorMaxwellian::Real3 temperature{0.0, 0.0, 0.0};
  ParticleGeneratorMaxwellian::Real3 pos{1.0, 2.0, 5.0};

  ParticleGeneratorMaxwellian gen{kind_idx, kind, mean_u, temperature};

  auto prt = gen.get(pos, {0.0, 0.0, 0.0});

  ASSERT_EQ(prt.kind, kind_idx);
  ASSERT_EQ(prt.w, w);
  ASSERT_EQ(prt.tag, 0);
  ASSERT_EQ(prt.u, mean_u); // zero temperature => exact velocity
  ASSERT_EQ(prt.x, pos);
}

// ======================================================================

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;

Grid_t* setupGrid()
{
  auto domain = Grid_t::Domain{{1, 8, 2},          // n grid points
                               {10.0, 80.0, 20.0}, // physical lengths
                               {0, 0, 0},          // location of lower corner
                               {1, 1, 1}};         // n patches

  auto bc =
    // FIXME wrong BCs
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 1.0, "i"};

  // --- generic setup
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 2;

  double dt = .1;
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

// TODO this template is dumb, just make a proper ctor
template <int PRT_COUNT, int KIND_IDX>
struct ParticleGenerator
{
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  psc::particle::Inject get(Real3 min_pos, Real3 pos_range)
  {
    Real uy = 2.0;
    if (PRT_COUNT > 0 && n_injected++ >= PRT_COUNT) {
      // Setting uy=0 means the particle won't enter the domain, and thus won't
      // be injected
      uy = 0.0;
    }

    Real3 x = min_pos + pos_range * Real3{0, .999, 0.};
    Real3 u{0.0, uy, 0.0};
    Real w = 1.0;
    psc::particle::Tag tag = 0;

    return {x, u, w, KIND_IDX, tag};
  }

  int n_injected = 0;
};

TEST(BoundaryInjectorTest, Integration1Particle)
{
  // ----------------------------------------------------------------------
  // setup

  PscParams psc_params;

  psc_params.nmax = 1;
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

  auto inject_particles =
    BoundaryInjector<ParticleGenerator<1, 1>, PscConfig::PushParticles>{{},
                                                                        grid};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks, marder,
                                          diagnostics, inject_particles);

  // ----------------------------------------------------------------------
  // set up initial conditions

  ASSERT_EQ(grid.n_patches(), 1);
  int p = 0;

  // ----------------------------------------------------------------------
  // run the simulation

  auto accessor = mprts.accessor();
  auto prts = accessor[p];

  ASSERT_EQ(prts.size(), 0);

  for (; grid.timestep_ < psc_params.nmax; grid.timestep_++) {
    psc.step();

    EXPECT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  ASSERT_EQ(prts.size(), 1);
}

TEST(BoundaryInjectorTest, IntegrationManyParticles)
{
  // ----------------------------------------------------------------------
  // setup

  PscParams psc_params;

  psc_params.nmax = 1;
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

  auto inject_particles =
    BoundaryInjector<ParticleGenerator<-1, 1>, PscConfig::PushParticles>{{},
                                                                         grid};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks, marder,
                                          diagnostics, inject_particles);

  // ----------------------------------------------------------------------
  // set up initial conditions

  ASSERT_EQ(grid.n_patches(), 1);
  int p = 0;

  // ----------------------------------------------------------------------
  // run the simulation

  auto accessor = mprts.accessor();
  auto prts = accessor[p];

  ASSERT_EQ(prts.size(), 0);

  for (; grid.timestep_ < psc_params.nmax; grid.timestep_++) {
    psc.step();

    EXPECT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  ASSERT_GT(prts.size(), 1);
}

TEST(BoundaryInjectorTest, IntegrationManySpecies)
{
  // ----------------------------------------------------------------------
  // setup

  PscParams psc_params;

  psc_params.nmax = 1;
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

  auto inject_electrons =
    BoundaryInjector<ParticleGenerator<-1, 0>, PscConfig::PushParticles>{{},
                                                                         grid};
  auto inject_ions =
    BoundaryInjector<ParticleGenerator<-1, 1>, PscConfig::PushParticles>{{},
                                                                         grid};

  auto inject = make_composite(inject_electrons, inject_ions);

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks, marder,
                                          diagnostics, inject);

  // ----------------------------------------------------------------------
  // set up initial conditions

  ASSERT_EQ(grid.n_patches(), 1);
  int p = 0;

  // ----------------------------------------------------------------------
  // run the simulation

  auto accessor = mprts.accessor();
  auto prts = accessor[p];

  ASSERT_EQ(prts.size(), 0);

  for (; grid.timestep_ < psc_params.nmax; grid.timestep_++) {
    psc.step();

    EXPECT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  bool found_electrons = false;
  bool found_ions = false;
  for (auto prt : prts) {
    found_electrons |= prt.kind() == 0;
    found_ions |= prt.kind() == 1;
  }

  ASSERT_TRUE(found_electrons);
  ASSERT_TRUE(found_ions);
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
