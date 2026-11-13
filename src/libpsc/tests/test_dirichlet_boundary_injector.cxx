#include <gtest/gtest.h>

#include "test_common.hxx"
#include "open_boundary_test_grid.hxx"

#include "psc.hxx"
#include "../psc_config.hxx"
#include "../psc_particle_injectors/dirichlet_boundary_injector.hxx"

TEST(DirichletBoundaryInjectorTest, ParticleGeneratorMaxwellianTest)
{
  int kind_idx = 15;
  Grid_t::Kind kind{1.0, 1836.0, "ion"};
  ParticleGeneratorMaxwellian::Real w = 1.0;
  ParticleGeneratorMaxwellian::Real3 mean_v{0.0, .05, .15};
  ParticleGeneratorMaxwellian::Real3 temperature{0.0, 0.0, 0.0};
  ParticleGeneratorMaxwellian::Real3 pos{1.0, 2.0, 5.0};

  ParticleGeneratorMaxwellian gen{kind_idx, kind, mean_v, temperature};

  auto prt = gen.get(pos, {0.0, 0.0, 0.0});

  ASSERT_EQ(prt.kind, kind_idx);
  ASSERT_EQ(prt.w, w);
  ASSERT_EQ(prt.tag, 0);
  ASSERT_EQ(prt.u / sqrt(1 + prt.u.mag2()),
            mean_v); // zero temperature => exact velocity
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
using psc::Axis;

// 8 cells along `d` (y or z), open there; 2 cells along the other axis,
// periodic
Grid_t* setupGrid(Axis d)
{
  return makeOpenBoundaryTestGrid(
    Int3{1, 2, 2}.with_component(d, 8),
    Vec3<bool>{false, false, false}.with_component(d, true), 1.);
}

struct ParticleGenerator
{
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  /// places particles at fraction `x_frac` of the given range, with momentum
  /// `u`; by default just inside the lower y boundary, moving in
  ParticleGenerator(int max_n_injected, int kind_idx,
                    Real3 x_frac = {0, .999, 0}, Real3 u = {0, 2, 0})
    : max_n_injected(max_n_injected), kind_idx(kind_idx), x_frac{x_frac}, u{u}
  {}

  psc::particle::Inject get(Real3 min_pos, Real3 pos_range)
  {
    Real3 u = this->u;
    if (max_n_injected > 0 && n_injected++ >= max_n_injected) {
      // Setting u=0 means the particle won't enter the domain, and thus won't
      // be injected
      u = {0, 0, 0};
    }

    Real3 x = min_pos + pos_range * x_frac;
    Real w = 1.0;
    psc::particle::Tag tag = 0;

    return {x, u, w, kind_idx, tag};
  }

  int n_injected = 0;
  int max_n_injected;
  int kind_idx;
  Real3 x_frac;
  Real3 u;
};

// generator for one particle entering through the given boundary normal to `d`
ParticleGenerator enteringGenerator(Axis d, LoHi lohi)
{
  bool lo = lohi == LoHi::Lo;
  return {1, KIND_ION,
          ParticleGenerator::Real3{0, 0, 0}.with_component(d, lo ? .999 : .001),
          ParticleGenerator::Real3{0, 0, 0}.with_component(d, lo ? 2 : -2)};
}

// Injects one particle through the given boundary and runs two steps, checking
// continuity and Gauss's law every step.
void integration_1_particle(Axis d, LoHi lohi)
{
  PscParams psc_params;

  psc_params.nmax = 2;
  psc_params.stats_every = 1;
  psc_params.cfl = .75;

  auto grid_ptr = setupGrid(d);
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
    new DirichletBoundaryInjector<ParticleGenerator, PscConfig::PushParticles>(
      enteringGenerator(d, lohi), d, lohi));

  ASSERT_EQ(grid.n_patches(), 1);
  ASSERT_EQ(mprts.size(), 0);

  psc.pre_first_step();
  for (; grid.timestep_ < psc_params.nmax;) {
    psc.step();

    ASSERT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    ASSERT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  ASSERT_EQ(mprts.size(), 1);
}

TEST(DirichletBoundaryInjectorTest, Integration1Particle)
{
  integration_1_particle(Axis::Y, LoHi::Lo);
}

TEST(DirichletBoundaryInjectorTest, IntegrationManyParticles)
{
  // ----------------------------------------------------------------------
  // setup

  PscParams psc_params;

  psc_params.nmax = 2;
  psc_params.stats_every = 1;
  psc_params.cfl = .75;

  auto grid_ptr = setupGrid(Axis::Y);
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
    new DirichletBoundaryInjector<ParticleGenerator, PscConfig::PushParticles>(
      ParticleGenerator(-1, 1), Axis::Y, LoHi::Lo));

  // ----------------------------------------------------------------------
  // set up initial conditions

  ASSERT_EQ(grid.n_patches(), 1);
  int p = 0;

  // ----------------------------------------------------------------------
  // run the simulation

  auto accessor = mprts.accessor();
  auto prts = accessor[p];

  ASSERT_EQ(prts.size(), 0);

  psc.pre_first_step();
  for (; grid.timestep_ < psc_params.nmax;) {
    psc.step();

    ASSERT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    ASSERT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  ASSERT_GT(prts.size(), 1);
}

TEST(DirichletBoundaryInjectorTest, IntegrationManySpecies)
{
  // ----------------------------------------------------------------------
  // setup

  PscParams psc_params;

  psc_params.nmax = 2;
  psc_params.stats_every = 1;
  psc_params.cfl = .75;

  auto grid_ptr = setupGrid(Axis::Y);
  auto& grid = *grid_ptr;

  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  ChecksParams checks_params{};
  checks_params.continuity.check_interval = 1;
  checks_params.gauss.check_interval = 1;
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  using Injector =
    DirichletBoundaryInjector<ParticleGenerator, PscConfig::PushParticles>;
  auto inject_electrons = Injector{ParticleGenerator(-1, 0), Axis::Y, LoHi::Lo};
  auto inject_ions = Injector{ParticleGenerator(-1, 1), Axis::Y, LoHi::Lo};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);

  psc.add_injector(&inject_ions);
  psc.add_injector(&inject_electrons);

  // ----------------------------------------------------------------------
  // set up initial conditions

  ASSERT_EQ(grid.n_patches(), 1);
  int p = 0;

  // ----------------------------------------------------------------------
  // run the simulation

  auto accessor = mprts.accessor();
  auto prts = accessor[p];

  ASSERT_EQ(prts.size(), 0);

  psc.pre_first_step();
  for (; grid.timestep_ < psc_params.nmax;) {
    psc.step();

    ASSERT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    ASSERT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
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

// Preacceleration sampling the normal E at the lower boundary must average the
// ghost and edge values, which straddle it.
void preaccelerate_samples_boundary(Axis d)
{
  auto grid_ptr = setupGrid(d);
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  int other = d == Axis::Y ? 2 : 1;
  auto F = make_Fields3d<Dim>(mflds[0]);
  for (int t = -1; t < grid.ldims[other] + 1; t++) {
    Int3 i = Int3{0, 0, 0}.with_component(other, t);
    i[d] = -1;
    F(EX + d, i[0], i[1], i[2]) = .2;
    i[d] = 0;
    F(EX + d, i[0], i[1], i[2]) = .1;
  }

  // the generator's second particle (u = 0) is kicked in too
  DirichletBoundaryInjector<ParticleGenerator, PscConfig::PushParticles>
    injector{enteringGenerator(d, LoHi::Lo), d, LoHi::Lo};
  injector.preaccelerate_method = PreaccelerateMethod::NormalE;
  injector.preaccelerate_time = 1.; // q/m = 1, so the kick is du = E
  injector.preaccelerate_sample_boundary = true;
  injector.inject(mprts, mflds);

  auto accessor = mprts.accessor();
  auto prts = accessor[0];
  ASSERT_EQ(prts.size(), 2);
  std::vector<double> u{prts[0].u()[d], prts[1].u()[d]};
  std::sort(u.begin(), u.end());
  EXPECT_NEAR(u[0], .15, 1e-6);
  EXPECT_NEAR(u[1], 2. + .15, 1e-6);
}

TEST(DirichletBoundaryInjectorTest, PreaccelerateSamplesBoundary)
{
  preaccelerate_samples_boundary(Axis::Y);
}

// The same at z boundaries.

TEST(DirichletBoundaryInjectorZTest, Integration1ParticleLo)
{
  integration_1_particle(Axis::Z, LoHi::Lo);
}

TEST(DirichletBoundaryInjectorZTest, Integration1ParticleHi)
{
  integration_1_particle(Axis::Z, LoHi::Hi);
}

TEST(DirichletBoundaryInjectorZTest, UpperDepositStartsAtOuterGhostFace)
{
  // The checks don't cover the upper boundary node, so check the current
  // directly: depositing from the ghost cell's outer face (z = ldims + 1), all
  // of the charge crosses the face at ldims + 1/2, and the fraction that ends
  // up on node ldims - 1 crosses the face at ldims - 1/2.
  auto grid_ptr = setupGrid(Axis::Z);
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  DirichletBoundaryInjector<ParticleGenerator, PscConfig::PushParticles>
    injector{enteringGenerator(Axis::Z, LoHi::Hi), Axis::Z, LoHi::Hi};
  injector.inject(mprts, mflds);

  auto accessor = mprts.accessor();
  auto prts = accessor[0];
  ASSERT_EQ(prts.size(), 1);
  int nz = grid.ldims[2];
  double z_final = prts[0].x()[2];
  ASSERT_GT(z_final, nz - 1);

  // J_z at index k is on the edge between nodes k and k + 1
  auto F = make_Fields3d<Dim>(mflds[0]);
  double jz_outer = 0, jz_inner = 0;
  for (int j = -1; j < grid.ldims[1] + 1; j++) {
    jz_outer += F(JXI + 2, 0, j, nz);
    jz_inner += F(JXI + 2, 0, j, nz - 1);
  }
  ASSERT_NE(jz_outer, 0.);
  EXPECT_NEAR(jz_inner / jz_outer, nz - z_final, 1e-10);
}

TEST(DirichletBoundaryInjectorZTest, PreaccelerateSamplesBoundary)
{
  preaccelerate_samples_boundary(Axis::Z);
}

TEST(DirichletBoundaryInjectorZTest, NormalPushIsRelativistic)
{
  // the particle moves only along the normal, by dt * u_z / gamma
  auto grid_ptr = setupGrid(Axis::Z);
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  ParticleGenerator generator{1, KIND_ION, {0, .5, .999}, {0, .5, 2}};
  DirichletBoundaryInjector<ParticleGenerator, PscConfig::PushParticles>
    injector{generator, Axis::Z, LoHi::Lo};
  injector.inject(mprts, mflds);

  auto accessor = mprts.accessor();
  auto prts = accessor[0];
  ASSERT_EQ(prts.size(), 1);
  double gamma = std::sqrt(1 + sqr(.5) + sqr(2.));
  EXPECT_NEAR(prts[0].x()[1], .5, 1e-12);
  EXPECT_NEAR(prts[0].x()[2], -.001 + grid.dt * 2. / gamma, 1e-12);
}

TEST(DirichletBoundaryInjectorTest, SinglePrecision)
{
  using PscConfigSingle = PscConfig1vbecSingle<Dim>;
  auto grid_ptr = setupGrid(Axis::Y);
  auto& grid = *grid_ptr;
  PscConfigSingle::MfieldsState mflds{grid};
  PscConfigSingle::Mparticles mprts{grid};

  DirichletBoundaryInjector<ParticleGenerator, PscConfigSingle::PushParticles>
    injector{enteringGenerator(Axis::Y, LoHi::Lo), Axis::Y, LoHi::Lo};
  injector.inject(mprts, mflds);

  EXPECT_EQ(mprts.size(), 1);
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
