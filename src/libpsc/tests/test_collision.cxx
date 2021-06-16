
#include "gtest/gtest.h"

#include "testing.hxx"
#include "../libpsc/psc_collision/psc_collision_impl.hxx"
#include "psc_particles_single.h"
#include "psc_fields_single.h"

#ifdef USE_CUDA
#include "../libpsc/cuda/collision_cuda_impl.hxx"
#endif

struct ParticleTest
{
  using real_t = double;
  using Real3 = Vec3<real_t>;

  Real3 u;
  real_t q;
  real_t m;
};

class MparticlesTest
{
public:
  using real_t = ParticleTest::real_t;

  real_t prt_q(const ParticleTest& prt) const { return prt.q; }
  real_t prt_m(const ParticleTest& prt) const { return prt.m; }
};

TEST(BinaryCollision, Test1)
{
  double eps = 1e-14;

  ParticleTest prt1{{1., 0., 0.}, 1., 1.};
  ParticleTest prt2{{0., 0., 0.}, 1., 1.};

  RngFake rng;
  MparticlesTest mprts;
  BinaryCollision<MparticlesTest, ParticleTest> bc(mprts);

  double nudt = bc(prt1, prt2, .1, rng);

  EXPECT_NEAR(prt1.u[0] + prt2.u[0], 1., eps);
  EXPECT_NEAR(prt1.u[1] + prt2.u[1], 0., eps);
  EXPECT_NEAR(prt1.u[2] + prt2.u[2], 0., eps);

#if 0
  printf("prt1: %g %g %g\n", prt1.u[0], prt1.u[1], prt1.u[2]);
  printf("prt2: %g %g %g\n", prt2.u[0], prt2.u[1], prt2.u[2]);
#endif
}

// ======================================================================
// make_psc
//
// FIXME, duplicated in various testing environments

template <typename dim>
static Grid_t& make_psc(const Grid_t::Kinds& kinds)
{
  Int3 gdims = {16, 16, 16};
  Int3 ibn = {2, 2, 2};
  Vec3<double> length = {160., 160., 160.};
  if (dim::InvarX::value) {
    gdims[0] = 1;
    ibn[0] = 0;
  }
  if (dim::InvarY::value) {
    gdims[1] = 1;
    ibn[1] = 0;
  }
  if (dim::InvarZ::value) {
    gdims[2] = 1;
    ibn[2] = 0;
  }

  auto grid_domain = Grid_t::Domain{gdims, length};
  auto grid_bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 200;
  auto coeff = Grid_t::Normalization{norm_params};
  return *new Grid_t{grid_domain, grid_bc, kinds, coeff, 1.};
}

// ======================================================================
// CollisionTest

template <typename DIM, typename COLLISION>
struct CollisionTestConfig
{
  using dim = DIM;
  using Collision = COLLISION;
  using Mparticles = typename Collision::Mparticles;
};

template <typename T>
class CollisionTest : public ::testing::Test
{};

using CollisionTestConfigSingle = CollisionTestConfig<
  dim_yz,
  CollisionHost<MparticlesSingle, MfieldsStateSingle, MfieldsSingle, RngFake>>;
using CollisionTestConfigDouble = CollisionTestConfig<
  dim_yz,
  CollisionHost<MparticlesDouble, MfieldsStateDouble, MfieldsC, RngFake>>;
#ifdef USE_CUDA
using CollisionTestConfigCuda =
  CollisionTestConfig<dim_yz,
                      CollisionCuda<MparticlesCuda<BS144>, RngStateFake>>;
#endif

using CollisionTestTypes =
  ::testing::Types<CollisionTestConfigSingle, CollisionTestConfigDouble
#ifdef xUSE_CUDA // FIXME
                   ,
                   CollisionTestConfigCuda
#endif
                   >;

TYPED_TEST_SUITE(CollisionTest, CollisionTestTypes);

TYPED_TEST(CollisionTest, Test1)
{
  using Config = TypeParam;
  using dim = typename Config::dim;
  using Mparticles = typename Config::Mparticles;
  using Collision = typename Config::Collision;
  const typename Mparticles::real_t eps = 1e-5;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  const auto& grid = make_psc<dim>(kinds);

  // init particles
  auto prt0 = psc::particle::Inject{{5., 5., 5.}, {1., 0., 0.}, 1., 0};
  auto prt1 = psc::particle::Inject{{5., 5., 5.}, {0., 0., 0.}, 1., 0};

  Mparticles mprts{grid};
  {
    auto inj = mprts.injector();
    auto injector = inj[0];
    injector(prt0);
    injector(prt1);
  }

  auto collision = Collision(grid, 1, 1.);

  collision(mprts);

  auto accessor = mprts.accessor();
  auto it = accessor[0].begin();
  auto prtf0 = *it++;
  auto prtf1 = *it++;
  EXPECT_NEAR(prtf0.u()[0] + prtf1.u()[0], 1., eps);
  EXPECT_NEAR(prtf0.u()[1] + prtf1.u()[1], 0., eps);
  EXPECT_NEAR(prtf0.u()[2] + prtf1.u()[2], 0., eps);

  // depends on random numbers, but for RngFake, we know
  EXPECT_NEAR(prtf0.u()[0], 0.96226911, eps);
  EXPECT_NEAR(prtf0.u()[1], 0., eps);
  EXPECT_NEAR(std::abs(prtf0.u()[2]), 0.17342988, eps);
  EXPECT_NEAR(prtf1.u()[0], 0.03773088, eps);
  EXPECT_NEAR(prtf1.u()[1], -0., eps);
  EXPECT_NEAR(std::abs(prtf1.u()[2]), 0.17342988, eps);
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
