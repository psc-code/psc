
#include "grid.hxx"
#include "fields.hxx"
#include "psc_particles_double.h"
#include "psc_fields_c.h"
#include "push_particles.hxx"
#include "../libpsc/psc_push_particles/push_config.hxx"
#include "../libpsc/psc_push_particles/push_dispatch.hxx"
#include "setup_fields.hxx"

#include "gtest/gtest.h"

// Rng hackiness

#include "../vpic/PscRng.h"

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

// ======================================================================
// class PushParticlesTest

struct PushParticlesTest : ::testing::Test
{
  using Mfields = MfieldsC;
  using Mparticles = MparticlesDouble;
  using PushParticles = PushParticles__<Config2nd<dim_yz>>;
  
  std::unique_ptr<Grid_t> grid_;

  const double L = 1e10;

  void SetUp()
  {
    auto domain = Grid_t::Domain{{2, 2, 2}, {L, L, L}};
    grid_.reset(new Grid_t{domain});
  }

};

// ======================================================================
// Accel test

TEST_F(PushParticlesTest, Accel)
{
  const int n_prts = 131;
  const int n_steps = 10;
  const Mparticles::real_t eps = 1e-6;

  // init fields
  auto mflds = Mfields{*grid_, NR_FIELDS, {2, 2, 2}};
  SetupFields<Mfields>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  // init particles
  grid_->kinds.push_back(Grid_t::Kind(1., 1., "test_species"));

  RngPool rngpool;
  Rng *rng = rngpool[0];
  auto n_prts_by_patch = std::vector<uint>{n_prts};

  auto mprts = Mparticles{*grid_};
  mprts.reserve_all(n_prts_by_patch.data());
  mprts.resize_all(n_prts_by_patch.data());
  for (auto& prt : mprts[0]) {
    prt.xi = rng->uniform(0, L);
    prt.yi = rng->uniform(0, L);
    prt.zi = rng->uniform(0, L);
    prt.qni_wni_ = 1.;
    prt.pxi = 0.;
    prt.pyi = 0.;
    prt.pzi = 0.;
    prt.kind_ = 0;
  }

  // run test
  PushParticles pushp_;
  for (int n = 0; n < n_steps; n++) {
    pushp_.push_mprts(mprts, mflds);

    for (auto& prt : mprts[0]) {
      EXPECT_NEAR(prt.pxi, 1*(n+1), eps);
      EXPECT_NEAR(prt.pyi, 2*(n+1), eps);
      EXPECT_NEAR(prt.pzi, 3*(n+1), eps);
    }
  }
}

// ======================================================================
// Cyclo test

TEST_F(PushParticlesTest, Cyclo)
{
  const int n_prts = 131;
  const int n_steps = 64;
  // the errors here are (substantial) truncation error, not
  // finite precision, and they add up
  // (but that's okay, if a reminder that the 6th order Boris would
  //  be good)
  const Mparticles::real_t eps = 1e-2;

  // init fields
  auto mflds = Mfields{*grid_, NR_FIELDS, {2, 2, 2}};
  SetupFields<Mfields>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case HZ: return 2. * M_PI / n_steps;
      default: return 0.;
      }
    });

  // init particles
  grid_->kinds.push_back(Grid_t::Kind(2., 1., "test_species"));

  RngPool rngpool;
  Rng *rng = rngpool[0];
  auto n_prts_by_patch = std::vector<uint>{n_prts};

  auto mprts = Mparticles{*grid_};
  mprts.reserve_all(n_prts_by_patch.data());
  mprts.resize_all(n_prts_by_patch.data());
  for (auto& prt : mprts[0]) {
    prt.xi = rng->uniform(0, L);
    prt.yi = rng->uniform(0, L);
    prt.zi = rng->uniform(0, L);
    prt.pxi = 1.; // gamma = 2
    prt.pyi = 1.;
    prt.pzi = 1.;
    prt.qni_wni_ = rng->uniform(0, 1.);;
    prt.kind_ = 0;
  }

  // run test
  PushParticles pushp_;
  for (int n = 0; n < n_steps; n++) {
    pushp_.push_mprts(mprts, mflds);

    double ux = (cos(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 cos(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uy = (sin(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 sin(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uz = 1.;
    for (auto& prt : mprts[0]) {
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
