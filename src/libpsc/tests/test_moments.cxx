
#include "gtest/gtest.h"

#include "testing.hxx"

using PushParticlesTestTypes = ::testing::Types<TestConfig1vbec3dSingleYZ
						,TestConfig1vbec3dSingle
#ifdef USE_CUDA
						,TestConfig1vbec3dCudaYZ
						,TestConfig1vbec3dCuda444
#endif
						>;

// FIXME, obviously this name is bad...
TYPED_TEST_CASE(PushParticlesTest, PushParticlesTestTypes);

TYPED_TEST(PushParticlesTest, Moment1)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using Moment_n = typename TypeParam::Moment_n;
  using Bnd = typename TypeParam::Bnd;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename Mfields::real_t;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts[0].injector();
    injector({{5., 5., 5.}, {0., 0., 1.}, 1., 0});
  }
  Moment_n moment_n{grid, grid.comm()};
  moment_n.run(mprts);
  auto& mres = moment_n.result();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	real_t val = mres[p](0, i,j,k);
	if (i == 0 && j == 0 & k == 0) {
	  EXPECT_NEAR(val, .005, eps) << "ijk " << i << " " << j << " " << k;
	} else {
	  EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
	}
	// if (val) {
	//   printf("ijk %d %d %d val %g\n", i, j, k, val);
	// }
      });
  }
}

TYPED_TEST(PushParticlesTest, Moment2) // FIXME, mostly copied
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using Moment_n = typename TypeParam::Moment_n;
  using Bnd = typename TypeParam::Bnd;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename Mfields::real_t;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts[0].injector();
    injector({{25., 5., 5.}, {0., 0., 1.}, 1., 0});
  }

  int i0 = 2;
  if (PushParticlesTest<TypeParam>::dim::InvarX::value) i0 = 0;
  
  Moment_n moment_n{grid, grid.comm()};
  moment_n.run(mprts);
  auto& mres = moment_n.result();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	real_t val = mres[p](0, i,j,k);
	if (i == i0 && j == 0 & k == 0) {
	  EXPECT_NEAR(val, .005, eps) << "ijk " << i << " " << j << " " << k;
	} else {
	  EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
	}
	// if (val) {
	//   printf("ijk %d %d %d val %g\n", i, j, k, val);
	// }
      });
  }
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
