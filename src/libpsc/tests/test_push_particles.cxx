
#include "gtest/gtest.h"

#include "testing.hxx"

#include "grid.hxx"
#include "fields.hxx"
#include "psc_particles_double.h"
#include "psc_fields_c.h"
#include "push_particles.hxx"
#include "setup_fields.hxx"
#include "setup_particles.hxx"

using PushParticlesTestTypes = ::testing::Types<TestConfig2ndDoubleYZ,
						TestConfig1vbec3dSingleYZ,
						TestConfig1vbec3dSingleXZ,
						//TestConfigVpic,
#ifdef USE_CUDA
						TestConfig1vbec3dCudaYZ,
						TestConfig1vbec3dCuda,
						TestConfig1vbec3dCuda444,
#endif
						TestConfig2ndDouble,
						TestConfig2ndSingle,
						TestConfig1vbec3dSingle>;

TYPED_TEST_CASE(PushParticlesTest, PushParticlesTestTypes);

// ======================================================================
// SingleParticle test

TYPED_TEST(PushParticlesTest, SingleParticle)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  const typename Mparticles::real_t eps = 1e-5;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init particle
  auto n_prts_by_patch = std::vector<uint>{1};

  Mparticles mprts{grid};
  {
    auto injector = mprts[0].injector();
    injector({{1., 0., 0.}, {}, 1., 0});
  }

  // check particle
  EXPECT_EQ(mprts[0].size(), 1);
  auto prt = *mprts[0].get().begin();
  EXPECT_NEAR(prt.position()[0], 1., eps);
  EXPECT_NEAR(prt.position()[1], 0., eps);
  EXPECT_NEAR(prt.position()[2], 0., eps);
  EXPECT_NEAR(prt.w(), 1., eps);
}

// ======================================================================
// vx, vy, vz

typename particle_inject::real_t vx(const particle_inject& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.u[0]) + sqr(prt.u[1]) + sqr(prt.u[2]));
  return gamma * prt.u[0];
}

typename particle_inject::real_t vy(const particle_inject& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.u[0]) + sqr(prt.u[1]) + sqr(prt.u[2]));
  return gamma * prt.u[1];
}

typename particle_inject::real_t vz(const particle_inject& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.u[0]) + sqr(prt.u[1]) + sqr(prt.u[2]));
  return gamma * prt.u[2];
}

// ======================================================================
// SingleParticlePushp1 test
//
// zero EM test, check position push given existing pzi

TYPED_TEST(PushParticlesTest, SingleParticlePushp1)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{5., 5., 5.}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  prt1.x[2] += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp2 test
//
// EZ = 1, check velocity push

TYPED_TEST(PushParticlesTest, SingleParticlePushp2)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return 2.;
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{5., 5., 5.}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  prt1.u[2] = 3.;
  prt1.x[2] += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp3 test
//
// EZ = z, check interpolation in z direction
// 1vbec doesn't actually interpolate EZ in the z direction,
// but by choosing the midpoint, that works out
// (but this test isn't very comprehensive)

TYPED_TEST(PushParticlesTest, SingleParticlePushp3)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return crd[2];
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{5., 5., 5.}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  prt1.u[2] = 6.;
  prt1.x[2] += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp4 test
//
// EZ = y, check interpolation in y direction

TYPED_TEST(PushParticlesTest, SingleParticlePushp4)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return crd[1];
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{5., 4., 5.}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  if (!Base::dim::InvarY::value) { prt1.u[2] = 5.; }
  this->push_x(prt0, prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp5 test
//
// EZ = x, check interpolation in x direction

TYPED_TEST(PushParticlesTest, SingleParticlePushp5)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return crd[0];
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{3., 5., 5.}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  prt1.u[2] = 4.;
  if (Base::dim::InvarX::value) { prt1.u[2] = 1.; }
  prt1.x[2] += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp6 test
//
// simply updating xi (in xyz, but not yz)

TYPED_TEST(PushParticlesTest, SingleParticlePushp6)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{1., 2., 3.}, {1., 1., 1.}, 1., 0};
  auto prt1 = prt0;
  if (!Base::dim::InvarX::value) prt1.x[0] += vx(prt1);
  if (!Base::dim::InvarY::value) prt1.x[1] += vy(prt1);
  if (!Base::dim::InvarZ::value) prt1.x[2] += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp7 test
//
// test with EZ = z in different block (not just lower left)

TYPED_TEST(PushParticlesTest, SingleParticlePushp7)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return crd[2];
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{151., 152., 155.}, {1., 1., 1.}, 1., 0};
  auto prt1 = prt0;
  prt1.u[2] = 156;
  this->push_x(prt0, prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp8 test
//
// test current deposition in simple case

TYPED_TEST(PushParticlesTest, SingleParticlePushp8)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{10., 10., 10.}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JZI, {1, 1, 1}, this->fnqz / this->dz * (xi1[2] - prt0.x[2]) },
    };
  }
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp9 test
//
// test current deposition, z move only, but crossing bnd

TYPED_TEST(PushParticlesTest, SingleParticlePushp9)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{10., 10., 19.5}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JZI, {1, 1, 1}, this->fnqz / this->dz * (20. - prt0.x[2]) },
      { JZI, {1, 1, 2}, this->fnqz / this->dz * (xi1[2] - 20.) },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp10 test
//
// test current deposition, y move only, but crossing bnd

TYPED_TEST(PushParticlesTest, SingleParticlePushp10)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{10., 19.5, 10.}, {0., 1., 0.}, 1., 0};
  auto prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    auto fnqy = .05;
    auto dy = 10.;
    curr_ref = {
      { JYI, {1, 1, 1}, fnqy / dy * (20. - prt0.x[1]) },
      { JYI, {1, 2, 1}, fnqy / dy * (xi1[1] - 20.) },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp11 test
//
// test current deposition, x move only

TYPED_TEST(PushParticlesTest, SingleParticlePushp11)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{10., 10., 10.}, {1., 0., 0.}, 1., 0};
  auto prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JXI, {1, 1, 1}, this->fnqx / this->dx * (xi1[0] - prt0.x[0]) },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp12 test
//
// test current deposition, yz move, no bnd crossing

TYPED_TEST(PushParticlesTest, SingleParticlePushp12)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{10., 10., 10.}, {0., 1., 1.}, 1., 0};
  auto prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JYI, {1, 1, 1}, 0.00280342 },
      { JYI, {1, 1, 2}, 8.333333e-05 },
      { JZI, {1, 1, 1}, 0.00280342 },
      { JZI, {1, 2, 1}, 8.333333e-05 },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp13 test
//
// test current deposition, yz move, cross z

TYPED_TEST(PushParticlesTest, SingleParticlePushp13)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{10., 19.5, 10.}, {0., 1., 1.}, 1., 0};
  auto prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JYI, {1, 1, 1}, 0.00243749 },
      { JZI, {1, 1, 1}, 6.25e-5 },
      { JYI, {1, 2, 1}, 0.00036592 },
      { JZI, {1, 2, 1}, 0.00282275 },
      { JYI, {1, 1, 2}, 6.25e-5 },
      { JYI, {1, 2, 2}, 2.08e-5 },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp14 test
//
// test current deposition, yz move, cross y

TYPED_TEST(PushParticlesTest, SingleParticlePushp14)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{10., 10., 19.5}, {0., 1., 1.}, 1., 0};
  auto prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JZI, {1, 1, 1}, 0.00243749 },
      { JYI, {1, 1, 1}, 6.25e-5 },
      { JZI, {1, 2, 1}, 6.25e-5 },
      { JZI, {1, 1, 2}, 0.00036592 },
      { JYI, {1, 1, 2}, 0.00282275 },
      { JZI, {1, 2, 2}, 2.08e-5 },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp15 test
//
// simple mv in z direction, cross block bnd

TYPED_TEST(PushParticlesTest, SingleParticlePushp15)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{5., 5., 39.5}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  this->push_x(prt0, prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp16 test
//
// simple mv in z direction, cross patch bnd

TYPED_TEST(PushParticlesTest, SingleParticlePushp16)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  auto prt0 = particle_inject{{5., 5., 159.5}, {0., 0., 1.}, 1., 0};
  auto prt1 = prt0;
  this->push_x(prt0, prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
