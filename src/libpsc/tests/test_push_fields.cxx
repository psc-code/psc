
#include "gtest/gtest.h"

#include "testing.hxx"

template<typename T>
struct PushFieldsTest : PushParticlesTest<T>
{
};

using PushFieldsTestTypes = ::testing::Types<TestConfig1vbec3dSingleYZ,
#ifdef USE_CUDA
					     TestConfig1vbec3dCudaYZ,
					     TestConfig1vbec3dCuda444,
#endif
					     TestConfig1vbec3dSingle>;

TYPED_TEST_CASE(PushFieldsTest, PushFieldsTestTypes);

// ======================================================================
// Pushf1

TYPED_TEST(PushFieldsTest, Pushf1)
{
  using Mparticles = typename TypeParam::Mparticles;
  using MfieldsState = typename TypeParam::MfieldsState;
  using dim = typename TypeParam::dim;
  using PushFields = typename TypeParam::PushFields;

  const typename Mparticles::real_t eps = 1e-2;

  this->make_psc({});
  const auto& grid = this->grid();

  const double kz = 2. * M_PI / grid.domain.length[2];
  
  // init fields
  auto mflds = MfieldsState{grid};
  SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case EY: return sin(kz*crd[2]);
      default: return 0.;
      }
    });

  // run test
  PushFields pushf_;
  pushf_.push_H(mflds, 1., dim{});

  auto flds = mflds[0];
  grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      double z = grid.patches[0].z_cc(k);
      // printf("ijk %d:%d:%d %g %g dt %g\n", i,j,k, flds(HX, i,j,k),
      // 	     kz*cos(kz*z), ppsc->dt);
      EXPECT_NEAR(flds(HX, i,j,k), kz*cos(kz*z), eps);
    });
}

TYPED_TEST(PushFieldsTest, Pushf2)
{
  using Mparticles = typename TypeParam::Mparticles;
  using MfieldsState = typename TypeParam::MfieldsState;
  using dim = typename TypeParam::dim;
  using PushFields = typename TypeParam::PushFields;

  const typename Mparticles::real_t eps = 1e-2;

  this->make_psc({});
  const auto& grid = this->grid();

  const double ky = 2. * M_PI / grid.domain.length[1];
  
  // init fields
  auto mflds = MfieldsState{grid};
  SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case HX: return cos(ky*crd[1]);
      default: return 0.;
      }
    });

  // run test
  PushFields pushf_;
  pushf_.push_E(mflds, 1., dim{});

  auto flds = mflds[0];
  grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      double y = grid.patches[0].y_cc(j);
      // printf("ijk %d:%d:%d %g %g dt %g\n", i,j,k, flds(HX, i,j,k),
      // 	     kz*cos(kz*z), ppsc->dt);
      EXPECT_NEAR(flds(EZ, i,j,k), ky*sin(ky*y), eps);
    });
}


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
