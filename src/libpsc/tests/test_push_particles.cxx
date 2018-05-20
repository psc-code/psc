
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
  SetupParticles<Mparticles>::setup_particles(mprts, n_prts_by_patch, [&](int p, int n) -> typename Mparticles::particle_t {
      typename Mparticles::particle_t prt{};
      prt.xi = 1.;
      prt.yi = 0.;
      prt.zi = 0.;
      prt.qni_wni_ = 1.;
      return prt;
    });

  // check particle
  for (auto& prt : make_getter(mprts)[0]) {
    EXPECT_NEAR(prt.xi, 1., eps);
    EXPECT_NEAR(prt.yi, 0., eps);
    EXPECT_NEAR(prt.zi, 0., eps);
    EXPECT_NEAR(prt.qni_wni_, 1., eps);
  }
}

// ======================================================================
// vx, vy, vz

template<typename particle_t>
typename particle_t::real_t vx(const particle_t& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.pxi) + sqr(prt.pyi) + sqr(prt.pzi));
  return gamma * prt.pxi;
}

template<typename particle_t>
typename particle_t::real_t vy(const particle_t& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.pxi) + sqr(prt.pyi) + sqr(prt.pzi));
  return gamma * prt.pyi;
}

template<typename particle_t>
typename particle_t::real_t vz(const particle_t& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.pxi) + sqr(prt.pyi) + sqr(prt.pzi));
  return gamma * prt.pzi;
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

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.zi += vz(prt1);
  
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

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 3.;
  prt1.zi += vz(prt1);
  
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

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 6.;
  prt1.zi += vz(prt1);
  
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

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 4.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 5.;
  prt1.zi = 5.980580;
  
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

  particle_t prt0, prt1;

  prt0.xi = 3.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 4.;
  if (Base::dim::InvarX::value) { prt1.pzi = 1.; }
  prt1.zi += vz(prt1);
  
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

  particle_t prt0, prt1;

  prt0.xi = 1.; prt0.yi = 2.; prt0.zi = 3.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 1.; prt0.pyi = 1.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  if (!Base::dim::InvarX::value) prt1.xi += vx(prt1);
  prt1.yi += vy(prt1);
  prt1.zi += vz(prt1);
  
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

  particle_t prt0, prt1;

  prt0.xi = 151.; prt0.yi = 152.; prt0.zi = 155.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 1.; prt0.pyi = 1.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 156;
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

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 10.; prt0.zi = 10.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JZI, {1, 1, 1}, this->fnqz / this->dz * (xi1[2] - prt0.zi) },
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

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 10.; prt0.zi = 19.5;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JZI, {1, 1, 1}, this->fnqz / this->dz * (20. - prt0.zi) },
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

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 19.5; prt0.zi = 10.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 1.; prt0.pzi = 0.;
  prt0.kind_ = 0;

  prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    auto fnqy = .05;
    auto dy = 10.;
    curr_ref = {
      { JYI, {1, 1, 1}, fnqy / dy * (20. - prt0.yi) },
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

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 10.; prt0.zi = 10.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 1.; prt0.pyi = 0.; prt0.pzi = 0.;
  prt0.kind_ = 0;

  prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JXI, {1, 1, 1}, this->fnqx / this->dx * (xi1[0] - prt0.xi) },
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

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 10.; prt0.zi = 10.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 1.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
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

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 19.5; prt0.zi = 10.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 1.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
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

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 10.; prt0.zi = 19.5;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 1.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
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

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 39.5;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
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

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 159.5;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  this->push_x(prt0, prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// PushParticlesTest2 is the same, but won't test cuda (because checks doesn't work)

template<typename T>
struct PushParticlesTest2 : PushParticlesTest<T>
{};

using PushParticlesTest2Types = ::testing::Types<TestConfig2ndDouble,
						 TestConfig2ndSingle,
						 TestConfig1vbec3dSingle>;

TYPED_TEST_CASE(PushParticlesTest2, PushParticlesTest2Types);

// ======================================================================
// Accel test

TYPED_TEST(PushParticlesTest2, Accel)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  using Checks = typename TypeParam::Checks;

  const int n_prts = 131;
  const int n_steps = 10;
  const typename Mparticles::real_t eps = 1e-5;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init fields
  auto mflds = Mfields{grid, NR_FIELDS, this->ibn};
  SetupFields<Mfields>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  // init particles
  RngPool rngpool;
  Rng *rng = rngpool[0];
  auto n_prts_by_patch = std::vector<uint>{n_prts};

  Mparticles mprts{grid};
  SetupParticles<Mparticles>::setup_particles(mprts, n_prts_by_patch, [&](int p, int n) -> typename Mparticles::particle_t {
      typename Mparticles::particle_t prt{};
      prt.xi = rng->uniform(0, this->L);
      prt.yi = rng->uniform(0, this->L);
      prt.zi = rng->uniform(0, this->L);
      prt.qni_wni_ = 1.;
      prt.pxi = 0.;
      prt.pyi = 0.;
      prt.pzi = 0.;
      prt.kind_ = 0;
      return prt;
    });

  // run test
  PushParticles pushp_;
  ChecksParams checks_params{};
  checks_params.continuity_threshold = 1e-10;
  checks_params.continuity_verbose = false;
  Checks checks_{grid, MPI_COMM_WORLD, checks_params};
  for (int n = 0; n < n_steps; n++) {
    //checks_.continuity_before_particle_push(mprts);
    pushp_.push_mprts(mprts, mflds);
    //checks_.continuity_after_particle_push(mprts, mflds);

    for (auto& prt : make_getter(mprts)[0]) {
      EXPECT_NEAR(prt.pxi, 1*(n+1), eps);
      EXPECT_NEAR(prt.pyi, 2*(n+1), eps);
      EXPECT_NEAR(prt.pzi, 3*(n+1), eps);
      prt.xi = this->L/2;
      prt.yi = this->L/2;
      prt.zi = this->L/2;
    }
  }
}

// ======================================================================
// Cyclo test

TYPED_TEST(PushParticlesTest2, Cyclo)
{
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
