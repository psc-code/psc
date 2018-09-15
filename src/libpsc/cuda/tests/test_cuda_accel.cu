
#include "grid.hxx"
#include "fields.hxx"
#include "cuda_mfields.h"
#include "cuda_mparticles.h"
#include "cuda_push_particles.cuh"
#include "push_particles_cuda_impl.hxx"

#include "cuda_test.hxx"

#include "../vpic/PscRng.h"

#include <memory>

#include "gtest/gtest.h"

// Rng hackiness

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

// enum hackiness

enum { // FIXME, duplicated
#if 0
  JXI, JYI, JZI,
  EX , EY , EZ ,
  HX , HY , HZ ,
#endif
  N_FIELDS = 9,
};

// profile hackiness

#include "mrc_profile.h"

struct prof_globals prof_globals; // FIXME

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 0;
}

using CudaMparticles = cuda_mparticles<BS144>;

// ======================================================================
// class PushMprtsTest

struct PushMprtsTest : TestBase<CudaMparticles>, ::testing::Test
{
  std::unique_ptr<Grid_t> grid_;

  RngPool rngpool;
  
  const double L = 1e10;

  void SetUp()
  {
    auto domain = Grid_t::Domain{{1, 4, 4}, {L, L, L}}; // FIXME, grid size needs to be at least a block; could use an assert to make sure that's the case...
    auto bc = GridBc{};
    auto kinds = Grid_t::Kinds{};
    auto norm = Grid_t::Normalization{};
    double dt = 1.;
    grid_.reset(new Grid_t{domain, bc, kinds, norm, dt});
  }

  // FIXME, convenient interfaces like make_cmflds, make_cmprts
  // should be available generally
  template<typename S>
  std::unique_ptr<cuda_mfields> make_cmflds(S set)
  {
    auto cmflds = std::unique_ptr<cuda_mfields>(new cuda_mfields(*grid_, N_FIELDS, { 0, 2, 2 }));

    fields_single_t flds = cmflds->get_host_fields();
    Fields3d<fields_single_t> F(flds);

    auto ldims = grid_->ldims;

    // FIXME, initializes some ghosts too many, but that doesn't really hurt...
    for (int k = 0; k <= ldims[2]; k++) {
      for (int j = 0; j <= ldims[1]; j++) {
	F(EX, 0,j,k) = set(EX);
	F(EY, 0,j,k) = set(EY);
	F(EZ, 0,j,k) = set(EZ);
	F(HX, 0,j,k) = set(HX);
	F(HY, 0,j,k) = set(HY);
	F(HZ, 0,j,k) = set(HZ);
      }
    }
    
    cmflds->copy_to_device(0, flds, 0, N_FIELDS);
    cmflds->dump("accel.fld.json");
    flds.dtor();
  
    return cmflds;
  }

};

// ======================================================================
// Accel test

TEST_F(PushMprtsTest, Accel)
{
  const int n_prts = 131;
  const int n_steps = 10;
  const CudaMparticles::real_t eps = 1e-6;

  // init fields
  auto cmflds = make_cmflds([&] (int m) -> cuda_mfields::real_t {
      switch(m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  // init particles
  Rng *rng = rngpool[0];

  grid_->kinds.push_back(Grid_t::Kind(1., 1., "test_species"));
  std::unique_ptr<CudaMparticles> cmprts(make_cmprts(*grid_, n_prts, [&](int i) -> cuda_mparticles_prt {
	using Real3 = cuda_mparticles_prt::Real3;
	auto prt = cuda_mparticles_prt{Real3(Vec3<double>{rng->uniform(0, L), rng->uniform(0, L), rng->uniform(0, L)}), {}, 1., 0};
	return prt;
      }));
  
  // run test
  for (int n = 0; n < n_steps; n++) {
    CudaPushParticles_<CudaConfig1vbec3d<dim_yz, BS144>>::push_mprts(cmprts.get(), cmflds.get());

    for (auto prt: cmprts->get_particles(0)) {
      EXPECT_NEAR(prt.p[0], 1*(n+1), eps);
      EXPECT_NEAR(prt.p[1], 2*(n+1), eps);
      EXPECT_NEAR(prt.p[2], 3*(n+1), eps);
    }
  }
}

// ======================================================================
// Cyclo test

TEST_F(PushMprtsTest, Cyclo)
{
  const int n_prts = 131;
  const int n_steps = 64;
  // the errors here are (substantial) truncation error, not
  // finite precision, and they add up
  // (but that's okay, if a reminder that the 6th order Boris would
  //  be good)
  const CudaMparticles::real_t eps = 1e-2;

  // init fields
  auto cmflds = make_cmflds([&] (int m) -> cuda_mfields::real_t {
      switch(m) {
      case HZ: return 2. * M_PI / n_steps;
      default: return 0.;
      }
    });

  // init particles
  Rng *rng = rngpool[0];

  grid_->kinds.push_back(Grid_t::Kind(2., 1., "test_species"));
  std::unique_ptr<CudaMparticles> cmprts(make_cmprts(*grid_, n_prts, [&](int i) -> cuda_mparticles_prt {
	using real_t = cuda_mparticles_prt::real_t;
	using Real3 = cuda_mparticles_prt::Real3;
	return {Real3(Vec3<double>({rng->uniform(0, L), rng->uniform(0, L), rng->uniform(0, L)})),
	        Real3{1., 1., 1.}, // gamma = 2
	        real_t(rng->uniform(0, 1.)),
		0};
      }));

  // run test
  for (int n = 0; n < n_steps; n++) {
    CudaPushParticles_<CudaConfig1vbec3d<dim_yz, BS144>>::push_mprts(cmprts.get(), cmflds.get());

    double ux = (cos(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 cos(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uy = (sin(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 sin(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uz = 1.;
    for (auto prt: cmprts->get_particles(0)) {
      EXPECT_NEAR(prt.p[0], ux, eps);
      EXPECT_NEAR(prt.p[1], uy, eps);
      EXPECT_NEAR(prt.p[2], uz, eps);
    }
  }
}

// ======================================================================
// main

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
