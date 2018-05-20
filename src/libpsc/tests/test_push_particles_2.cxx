
#include "gtest/gtest.h"

#include "testing.hxx"

using PushParticlesTestTypes = ::testing::Types<TestConfig2ndDoubleYZ,
#ifdef USE_CUDA
						TestConfig1vbec3dCudaYZ,
#endif
#if 0
						TestConfig1vbec3dCuda,
						TestConfig1vbec3dCuda444,
#endif
						TestConfig1vbec3dSingleYZ>;
						//TestConfig1vbec3dSingle>;
						//TestConfig2ndDouble>;

TYPED_TEST_CASE(PushParticlesTest, PushParticlesTestTypes);

// ======================================================================
// Accel test

TYPED_TEST(PushParticlesTest, Accel)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  using BndParticles = typename TypeParam::BndParticles;
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
  BndParticles bndp_{ppsc->mrc_domain_, grid};
  ChecksParams checks_params{};
  checks_params.continuity_threshold = 1e-10;
  checks_params.continuity_verbose = false;
  Checks checks_{grid, MPI_COMM_WORLD, checks_params};
  for (int n = 0; n < n_steps; n++) {
    //checks_.continuity_before_particle_push(mprts);
    pushp_.push_mprts(mprts, mflds);
    //checks_.continuity_after_particle_push(mprts, mflds);
    bndp_(mprts);

    for (auto& prt : make_getter(mprts)[0]) {
      EXPECT_NEAR(prt.pxi, 1*(n+1), eps);
      EXPECT_NEAR(prt.pyi, 2*(n+1), eps);
      EXPECT_NEAR(prt.pzi, 3*(n+1), eps);
    }
  }
}

// ======================================================================
// Cyclo test

TYPED_TEST(PushParticlesTest, Cyclo)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  using BndParticles = typename TypeParam::BndParticles;
  using Checks = typename TypeParam::Checks;

  const int n_prts = 131;
  const int n_steps = 64;
  // the errors here are (substantial) truncation error, not
  // finite precision, and they add up
  // (but that's okay, if a reminder that the 6th order Boris would
  //  be good)
  const typename Mparticles::real_t eps = 1e-2;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(2., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init fields
  auto mflds = Mfields{grid, NR_FIELDS, this->ibn};
  SetupFields<Mfields>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case HZ: return 2. * M_PI / n_steps;
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
      prt.pxi = 1.; // gamma = 2
      prt.pyi = 1.;
      prt.pzi = 1.;
      prt.qni_wni_ = rng->uniform(0, 1.);;
      prt.kind_ = 0;
      return prt;
    });

  // run test
  PushParticles pushp_;
  BndParticles bndp_{ppsc->mrc_domain_, grid};
  ChecksParams checks_params{};
  checks_params.continuity_threshold = 1e-10;
  checks_params.continuity_verbose = false;
  Checks checks_{grid, MPI_COMM_WORLD, checks_params};
  for (int n = 0; n < n_steps; n++) {
    //checks_.continuity_before_particle_push(mprts);
    pushp_.push_mprts(mprts, mflds);
    //checks_.continuity_after_particle_push(mprts, mflds);
    bndp_(mprts);

    double ux = (cos(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 cos(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uy = (sin(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 sin(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uz = 1.;
    for (auto& prt : make_getter(mprts)[0]) {
      EXPECT_NEAR(prt.pxi, ux, eps);
      EXPECT_NEAR(prt.pyi, uy, eps);
      EXPECT_NEAR(prt.pzi, uz, eps);
    }
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
