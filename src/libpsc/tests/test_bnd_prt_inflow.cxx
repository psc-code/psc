
#include <gtest/gtest.h>

#include "test_common.hxx"

#include "inflow.hxx"
#include "psc_particles_double.h"
#include "psc_particles_single.h"
#include "particle_with_id.h"
#include "setup_particles.hxx"
#include "pushp.hxx"
#include "dim.hxx"
#ifdef USE_CUDA
#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "../libpsc/cuda/mparticles_cuda.inl"
#endif
#include "particles_simple.inl"
#include <kg/io.h>

template <typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using MparticlesInflowTestTypes =
  ::testing::Types<Config<MparticlesSingle>,
                   Config<MparticlesSingle, MakeTestGridYZ>
#ifdef USE_CUDA
                   ,
                   Config<MparticlesCuda<BS144>, MakeTestGridYZ1>,
                   Config<MparticlesCuda<BS144>, MakeTestGridYZ>
#endif
                   >;

TYPED_TEST_SUITE(MparticlesInflowTest, MparticlesInflowTestTypes);

TEST(TestSetupParticlesInflow, SetupSimple)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, .1};
  Mparticles mprts{grid};

  SetupParticles<Mparticles> setup_particles(grid);

  Double3 pos = {5., 5., 5.};
  Double3 u = Double3{0.0, 0.0, 0.0};
  int tag = 0;
  psc_particle_np np = {0, 1.0, [u]() { return u; }, tag};
  double wni = 1.0;
  auto prt = setup_particles.setupParticle(np, pos, wni);

  EXPECT_EQ(prt.x, pos);
  EXPECT_EQ(prt.u, u);
  EXPECT_EQ(prt.w, wni);
  EXPECT_EQ(prt.kind, 0);
  EXPECT_EQ(prt.tag, tag);
}

TEST(TestSetupParticlesInflow, SetupMaxwellian)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, .1};
  Mparticles mprts{grid};

  SetupParticles<Mparticles> setup_particles(grid);

  Double3 pos = {5., 5., 5.};
  Double3 u = {0.0, 0.5, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  int tag = 0;

  psc_particle_npt npt = {0, 1.0, u, T, tag};
  auto p = setup_particles.createMaxwellian(npt);
  psc_particle_np np = {0, 1.0, p, tag};
  double wni = 1.0;
  auto prt = setup_particles.setupParticle(np, pos, wni);

  EXPECT_EQ(prt.x, pos);
  EXPECT_EQ(prt.u, u);
  EXPECT_EQ(prt.w, wni);
  EXPECT_EQ(prt.kind, 0);
  EXPECT_EQ(prt.tag, tag);
}

class TestInjector
{
public:
  void operator()(psc::particle::Inject prt) { prts.push_back(prt); }

  std::vector<psc::particle::Inject> prts;
};

static double half() { return 0.5; }

TEST(TestSetupParticlesInflow, Advance)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 pos = {5., -2., 5.};
  Double3 u = {0.0, 0.5, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles, dim_y> inflow(grid, npt, *half);

  auto prt = inflow.get_advanced_prt(pos, 1.0);

  EXPECT_NEAR(prt.x[0], 5., 1e-5);
  EXPECT_NEAR(prt.x[1], 2.47214, 1e-5);
  EXPECT_NEAR(prt.x[2], 5., 1e-5);
}

TEST(TestSetupParticlesInflow, InjectIntoCell)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 u = {0.0, 0.8, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles, dim_y> inflow(grid, npt, *half);

  TestInjector injector;
  inflow.inject_into_boundary_cell(grid, injector, {0, 0, 0});

  EXPECT_EQ(injector.prts.size(), prm.nicell);

  for (auto prt : injector.prts) {
    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 1.246950, 1e-5);
    EXPECT_NEAR(prt.x[2], 5., 1e-5);
  }
}

TEST(TestSetupParticlesInflow, InjectIntoCellFilter)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  // too slow to enter domain
  Double3 u = {0.0, 0.5, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles, dim_y> inflow(grid, npt, *half);

  TestInjector injector;
  inflow.inject_into_boundary_cell(grid, injector, {0, 0, 0});

  EXPECT_EQ(injector.prts.size(), 0);
}

TEST(TestSetupParticlesInflow, InjectIntoPatch)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 u = {0.0, 0.8, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles, dim_y> inflow(grid, npt, *half);

  TestInjector injector;
  inflow.inject_into_boundary_patch(grid, injector, grid.patches[0]);

  EXPECT_EQ(injector.prts.size(),
            domain.ldims[0] * domain.ldims[2] * prm.nicell);

  for (int i = 0; i < injector.prts.size(); i++) {
    auto prt = injector.prts[i];

    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 1.246950, 1e-5);

    if (i < prm.nicell) {
      EXPECT_NEAR(prt.x[2], 5., 1e-5);
    } else {
      EXPECT_NEAR(prt.x[2], 15., 1e-5);
    }
  }
}

TEST(TestSetupParticlesInflow, InjectIntoBoundaryY)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 2, 2}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 u = {0.0, 0.8, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles, dim_y> inflow(grid, npt, *half);

  std::vector<TestInjector> injectors = {TestInjector(), TestInjector(),
                                         TestInjector(), TestInjector()};
  inflow.inject_into_boundary(grid, injectors);

  EXPECT_EQ(injectors[0].prts.size(),
            domain.ldims[0] * domain.ldims[2] * prm.nicell);
  for (auto prt : injectors[0].prts) {
    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 1.246950, 1e-5);
    EXPECT_NEAR(prt.x[2], 5., 1e-5);
  }

  EXPECT_EQ(injectors[1].prts.size(), 0);

  EXPECT_EQ(injectors[2].prts.size(),
            domain.ldims[0] * domain.ldims[2] * prm.nicell);
  for (auto prt : injectors[2].prts) {
    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 1.246950, 1e-5);
    EXPECT_NEAR(prt.x[2], 15., 1e-5);
  }

  EXPECT_EQ(injectors[3].prts.size(), 0);
}

TEST(TestSetupParticlesInflow, InjectIntoBoundaryZ)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 2, 2}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 u = {0.0, 0.0, 0.8};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles, dim_z> inflow(grid, npt, *half);

  std::vector<TestInjector> injectors = {TestInjector(), TestInjector(),
                                         TestInjector(), TestInjector()};
  inflow.inject_into_boundary(grid, injectors);

  EXPECT_EQ(injectors[0].prts.size(),
            domain.ldims[0] * domain.ldims[1] * prm.nicell);
  for (auto prt : injectors[0].prts) {
    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 5., 1e-5);
    EXPECT_NEAR(prt.x[2], 1.246950, 1e-5);
  }

  EXPECT_EQ(injectors[1].prts.size(),
            domain.ldims[0] * domain.ldims[1] * prm.nicell);
  for (auto prt : injectors[1].prts) {
    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 15., 1e-5);
    EXPECT_NEAR(prt.x[2], 1.246950, 1e-5);
  }

  EXPECT_EQ(injectors[2].prts.size(), 0);

  EXPECT_EQ(injectors[3].prts.size(), 0);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
