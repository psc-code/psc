
#include "gtest/gtest.h"

#include "testing.hxx"

using PushParticlesTestTypes = ::testing::Types<TestConfig2ndDoubleYZ
						,TestConfig1vbec3dSingle
#ifdef USE_CUDA
						,TestConfig1vbec3dCudaYZ
						//,TestConfig1vbec3dCuda444
#endif
						>;
//TestConfig1vbec3dSingle>;
//TestConfig2ndDouble>;

TYPED_TEST_CASE(PushParticlesTest, PushParticlesTestTypes);

// ======================================================================
// Accel test

TYPED_TEST(PushParticlesTest, Accel)
{
  using Mparticles = typename TypeParam::Mparticles;
  using MfieldsState = typename TypeParam::MfieldsState;
  using PushParticles = typename TypeParam::PushParticles;
  using BndParticles = typename TypeParam::BndParticles;
  using Bnd = typename TypeParam::Bnd;
  using Checks = typename TypeParam::Checks;

  const int n_prts = 131;
  const int n_steps = 10;
  const typename Mparticles::real_t eps = 1e-5;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init fields
  auto mflds = MfieldsState{grid};
  SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
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

  Mparticles mprts{grid};
  for (int p = 0; p < grid.n_patches(); p++) {
    auto injector = mprts[p].injector();
    for (int n = 0; n < n_prts; n++) {
      injector({{rng->uniform(0, this->L), rng->uniform(0, this->L), rng->uniform(0, this->L)},
	        {}, 1., 0});
    }
  }

  // run test
  PushParticles pushp_;
  BndParticles bndp_{grid};
  Bnd bnd_{grid, this->ibn};
  ChecksParams checks_params{};
  checks_params.continuity_threshold = 1e-10;
  checks_params.continuity_verbose = false;
  Checks checks_{grid, MPI_COMM_WORLD, checks_params};
  for (int n = 0; n < n_steps; n++) {
    checks_.continuity_before_particle_push(mprts);
    pushp_.push_mprts(mprts, mflds);
    bndp_(mprts);
    bnd_.add_ghosts(mflds, JXI, JXI+3);
    bnd_.fill_ghosts(mflds, JXI, JXI+3);
    checks_.continuity_after_particle_push(mprts, mflds);

    for (auto prt : mprts[0].get()) {
      EXPECT_NEAR(prt.u()[0], 1*(n+1), eps);
      EXPECT_NEAR(prt.u()[1], 2*(n+1), eps);
      EXPECT_NEAR(prt.u()[2], 3*(n+1), eps);
    }
  }
}

// ======================================================================
// Cyclo test

TYPED_TEST(PushParticlesTest, Cyclo)
{
  using Mparticles = typename TypeParam::Mparticles;
  using MfieldsState = typename TypeParam::MfieldsState;
  using PushParticles = typename TypeParam::PushParticles;
  using BndParticles = typename TypeParam::BndParticles;
  using Bnd = typename TypeParam::Bnd;
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
  auto mflds = MfieldsState{grid};
  SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case HZ: return 2. * M_PI / n_steps;
      default: return 0.;
      }
    });

  // init particles
  RngPool rngpool;
  Rng *rng = rngpool[0];

  Mparticles mprts{grid};
  for (int p = 0; p < grid.n_patches(); p++) {
    auto injector = mprts[p].injector();
    for (int n = 0; n < n_prts; n++) {
      injector({{rng->uniform(0, this->L), rng->uniform(0, this->L), rng->uniform(0, this->L)},
	        {1., 1., 1.},
		rng->uniform(0., 1.), 0});
    }
  }

  // run test
  PushParticles pushp_;
  BndParticles bndp_{grid};
  Bnd bnd_{grid, this->ibn};
  ChecksParams checks_params{};
  checks_params.continuity_threshold = 1e-10;
  checks_params.continuity_verbose = false;
  Checks checks_{grid, MPI_COMM_WORLD, checks_params};
  for (int n = 0; n < n_steps; n++) {
    checks_.continuity_before_particle_push(mprts);
    pushp_.push_mprts(mprts, mflds);
    bndp_(mprts);
    bnd_.add_ghosts(mflds, JXI, JXI+3);
    bnd_.fill_ghosts(mflds, JXI, JXI+3);
    checks_.continuity_after_particle_push(mprts, mflds);

    double ux = (cos(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 cos(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uy = (sin(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 sin(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uz = 1.;
    for (auto prt : mprts[0].get()) {
      EXPECT_NEAR(prt.u()[0], ux, eps);
      EXPECT_NEAR(prt.u()[1], uy, eps);
      EXPECT_NEAR(prt.u()[2], uz, eps);
    }
  }
}

// ======================================================================
// main

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  pr_time_step_no_comm = prof_register("time step w/o comm", 1., 0, 0);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
