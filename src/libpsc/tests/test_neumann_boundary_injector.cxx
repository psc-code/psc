#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc.hxx"
#include "../psc_config.hxx"
#include "../psc_particle_injectors/neumann_boundary_injector.hxx"

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;

using psc::Axis;

// 8 cells along the boundary normal `d` (y or z), 2 along the other axis
Grid_t::Domain makeDomain(Axis d)
{
  Int3 gdims = Int3{1, 2, 2}.with_component(d, 8);
  return Grid_t::Domain{gdims,          // n grid points
                        Double3(gdims), // physical lengths
                        {0, 0, 0},      // location of lower corner
                        {1, 1, 1}};     // n patches
}

// open along `d`, periodic otherwise
Grid_t* setupGridWithDt(Axis d, double dt)
{
  auto domain = makeDomain(d);

  Int3 fld =
    Int3{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC}.with_component(
      d, BND_FLD_OPEN);
  Int3 prt =
    Int3{BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}.with_component(
      d, BND_PRT_OPEN);
  auto bc = psc::grid::BC{fld, fld, prt, prt};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 1.0, "i"};

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 1;

  Grid_t::Normalization norm{norm_params};

  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

Grid_t* setupGrid(Axis d, double cfl)
{
  return setupGridWithDt(d, cfl * courant_length(makeDomain(d)));
}

// For each (y, uy), places a stationary electron and a moving ion at the same
// position (to satisfy Gauss' law at t=0), runs the simulation for `nmax` steps
// with a NeumannBoundaryInjector at the given boundary, and returns the
// final particles' positions along `d`, sorted. Here y and uy are the position
// and momentum along `d`. Copies are made at `copy_depth`.
template <typename Resampler = NeumannResamplerNone>
std::vector<double> run(Axis d, LoHi lohi,
                        std::vector<std::pair<double, double>> ys_uys, int nmax,
                        Resampler resampler = {}, int copy_depth = 1)
{
  PscParams psc_params;
  psc_params.nmax = nmax;
  psc_params.stats_every = 1;
  psc_params.cfl = .75;

  auto grid_ptr = setupGrid(d, psc_params.cfl);
  auto& grid = *grid_ptr;

  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  ChecksParams checks_params{};
  checks_params.continuity.check_interval = 1;
  checks_params.gauss.check_interval = 1;
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);

  NeumannBoundaryInjector<PscConfig::PushParticles, Resampler> injector{
    d, lohi, resampler, copy_depth};
  psc.add_injector(&injector);

  {
    auto inj = mprts.injector()[0];
    for (auto [y, uy] : ys_uys) {
      auto x = Double3{0, .5, .5}.with_component(d, y);
      auto u = Double3{0, 0, 0}.with_component(d, uy);
      inj({x, {0, 0, 0}, 1, KIND_ELECTRON});
      inj({x, u, 1, KIND_ION});
    }
  }

  psc.pre_first_step();
  for (; grid.timestep_ < psc_params.nmax;) {
    psc.step();

    EXPECT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  std::vector<double> ys;
  auto accessor = mprts.accessor();
  for (auto prt : accessor[0]) {
    ys.push_back(prt.position()[d]);
  }
  std::sort(ys.begin(), ys.end());
  return ys;
}

// v = 2 / sqrt(5) ~ .894 and dt ~ .53, so particles move ~.47 cells per step

void inwards_lo(Axis d)
{
  auto ys = run(d, LoHi::Lo, {{.75, 2.}}, 1);
  ASSERT_EQ(ys.size(), 3);
  EXPECT_NEAR(ys[0], ys[2] - 1., 1e-10); // copy, one cell behind the ion
  EXPECT_EQ(ys[1], .75);                 // electron
  EXPECT_GT(ys[2], 1.);                  // ion
}

void inwards_hi(Axis d)
{
  auto ys = run(d, LoHi::Hi, {{7.25, -2.}}, 1);
  ASSERT_EQ(ys.size(), 3);
  EXPECT_LT(ys[0], 7.);                  // ion
  EXPECT_EQ(ys[1], 7.25);                // electron
  EXPECT_NEAR(ys[2], ys[0] + 1., 1e-10); // copy, one cell behind the ion
}

void depth_hi(Axis d)
{
  auto ys = run(d, LoHi::Hi, {{5.25, -2.}}, 1, {}, 3); // crosses 5: copied
  ASSERT_EQ(ys.size(), 3);
  EXPECT_LT(ys[0], 5.);                  // ion
  EXPECT_EQ(ys[1], 5.25);                // electron
  EXPECT_NEAR(ys[2], ys[0] + 3., 1e-10); // copy, three cells behind the ion
}

TEST(NeumannBoundaryInjectorTest, InwardsLo) { inwards_lo(Axis::Y); }

TEST(NeumannBoundaryInjectorTest, InwardsHi) { inwards_hi(Axis::Y); }

TEST(NeumannBoundaryInjectorTest, NoCopies)
{
  auto ys = run(Axis::Y, LoHi::Lo,
                {{.25, 1.},    // stays in edge cell
                 {.25, -2.},   // leaves domain
                 {7.25, -2.}}, // moves inwards at other boundary
                1);
  ASSERT_EQ(ys.size(), 5); // outgoing ion was dropped
}

TEST(NeumannBoundaryInjectorTest, ManySteps)
{
  // each copy is itself copied when it leaves the edge cell; resampling with a
  // high temperature checks that the copies still conserve charge
  std::unique_ptr<Grid_t> grid{setupGrid(Axis::Y, .75)};
  auto ys = run(Axis::Y, LoHi::Lo, {{.75, 2.}}, 6,
                NeumannResamplerMaxwellian{grid->kinds, {9., 9.}});
  ASSERT_GT(ys.size(), 3);
}

// Copies made deeper than the edge cell (copy_depth > 1).

TEST(NeumannBoundaryInjectorDepthTest, Lo)
{
  // staggered so that a copy of the first ion would sit elsewhere
  auto ys = run(Axis::Y, LoHi::Lo,
                {{.6, 2.},    // crosses y = 1: not copied
                 {2.75, 2.}}, // crosses y = 3: copied
                1, {}, 3);
  ASSERT_EQ(ys.size(), 5);
  EXPECT_LT(ys[0], 1.);                  // copy, three cells behind the ion
  EXPECT_EQ(ys[1], .6);                  // electron
  EXPECT_GT(ys[2], 1.);                  // first ion
  EXPECT_EQ(ys[3], 2.75);                // electron
  EXPECT_NEAR(ys[0], ys[4] - 3., 1e-10); // second ion
}

TEST(NeumannBoundaryInjectorDepthTest, Hi) { depth_hi(Axis::Y); }

TEST(NeumannBoundaryInjectorDepthTest, ManySteps)
{
  // copies of copies, resampled hot, still conserve charge
  std::unique_ptr<Grid_t> grid{setupGrid(Axis::Y, .75)};
  auto ys = run(Axis::Y, LoHi::Lo, {{2.75, 2.}}, 10,
                NeumannResamplerMaxwellian{grid->kinds, {9., 9.}}, 3);
  ASSERT_GT(ys.size(), 3);
}

// The same at z boundaries.

TEST(NeumannBoundaryInjectorZTest, InwardsLo) { inwards_lo(Axis::Z); }

TEST(NeumannBoundaryInjectorZTest, InwardsHi) { inwards_hi(Axis::Z); }

TEST(NeumannBoundaryInjectorZTest, DepthHi) { depth_hi(Axis::Z); }

TEST(NeumannBoundaryInjectorTest, CopyAcrossPeriodicBoundary)
{
  // the ion moves inward and across the periodic z boundary in the same step,
  // so its copy starts outside the patch in z, until particle exchange wraps
  // it; injecting it must not move it
  PscParams psc_params;
  psc_params.nmax = 1;

  auto grid_ptr = setupGridWithDt(Axis::Y, 1.);
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};
  ChecksParams checks_params{};
  checks_params.continuity.check_interval = 1;
  checks_params.gauss.check_interval = 1;
  Checks checks{grid, MPI_COMM_WORLD, checks_params};
  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);
  psc.add_injector(
    new NeumannBoundaryInjector<PscConfig::PushParticles>{Axis::Y, LoHi::Lo});

  {
    auto injector = mprts.injector();
    auto inj = injector[0];
    inj({{.5, .9, .1}, {0, 0, 0}, 1, KIND_ELECTRON});
    inj({{.5, .9, .1}, {0, 1, -1}, 1, KIND_ION});
  }

  psc.pre_first_step();
  psc.step();

  ASSERT_EQ(mprts.size(), 3);
  EXPECT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
  EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);

  // the copy is one cell behind the ion
  auto accessor = mprts.accessor();
  std::vector<Vec3<double>> ions;
  for (auto prt : accessor[0]) {
    if (prt.kind() == KIND_ION) {
      ions.push_back(Vec3<double>(prt.x()));
    }
  }
  ASSERT_EQ(ions.size(), 2);
  EXPECT_NEAR(std::abs(ions[0][1] - ions[1][1]), 1., 1e-10);
  EXPECT_NEAR(ions[0][2], ions[1][2], 1e-10);
}

TEST(NeumannBoundaryInjectorTest, SinglePrecision)
{
  // an ion just pushed out of the edge cell, at y = 1.2 from y ~ .73
  using PscConfigSingle = PscConfig1vbecSingle<Dim>;
  auto grid_ptr = setupGrid(Axis::Y, .75);
  auto& grid = *grid_ptr;
  PscConfigSingle::MfieldsState mflds{grid};
  PscConfigSingle::Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{.5, 1.2, .5}, {0, 2, 0}, 1, KIND_ION});
  }

  NeumannBoundaryInjector<PscConfigSingle::PushParticles> injector{Axis::Y,
                                                                   LoHi::Lo};
  injector.inject(mprts, mflds);

  EXPECT_EQ(mprts.size(), 2);
}

TEST(NeumannResamplerTest, Rotate)
{
  NeumannResamplerRotate resampler;
  psc::particle::Inject prt{{0, 0, 0}, {.3, .5, .4}, 1, KIND_ION};
  auto u = resampler.resample(prt, 1);
  EXPECT_EQ(u[1], .5);
  EXPECT_NEAR(sqr(u[0]) + sqr(u[2]), .25, 1e-12);
  EXPECT_NE(u[0], .3);
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
