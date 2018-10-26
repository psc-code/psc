
#include <gtest/gtest.h>

#include "psc_config.h"
#include "test_common.hxx"
#include "mpi.h"

// FIXME, the general tests should be moved -> test_mparticles,
// and the real cuda ones to test_cuda_mparticle?

#ifdef USE_CUDA
#include "../libpsc/cuda/psc_particles_cuda.h"

template<typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using MparticlesCudaTestTypes = ::testing::Types<Config<MparticlesCuda<BS144>, MakeTestGridYZ>>;

TYPED_TEST_CASE(MparticlesCudaTest, MparticlesCudaTestTypes);

// ======================================================================
// MparticlesCudaTest

template<typename T>
struct MparticlesCudaTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using MakeGrid = typename T::MakeGrid;

  MparticlesCudaTest()
    : grid_{MakeGrid{}()}
  {
    grid_.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  }

  Mparticles mk_mprts()
  {
    Mparticles mprts(grid_);
    mprts.define_species("test_species", 1., 1., 100, 10,
			 10, 0);
    return mprts;
  }

private:
  Grid_t grid_;
};

// -----------------------------------------------------------------------
// Constructor

TYPED_TEST(MparticlesCudaTest, Constructor)
{
  auto mprts = this->mk_mprts();
}

// -----------------------------------------------------------------------
// Inject

TYPED_TEST(MparticlesCudaTest, Inject)
{
  using Mparticles = typename TypeParam::Mparticles;
  const int n_prts = 1;
  
  auto mprts = this->mk_mprts();

  EXPECT_EQ(mprts.n_patches(), 4);

  int nn = 0;
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto injector = mprts[p].injector();
    auto& patch = mprts.grid().patches[p];
    for (int n = 0; n < n_prts; n++) {
      particle_inject prt = {};
      auto x = .5 * (patch.xb + patch.xe);
      int kind = 0;
      // use weight to store particle number for testing
      injector({{x[0], x[1], x[2]}, {}, double(nn), kind});
      nn++;
    }
  }

  EXPECT_EQ(mprts.get_n_prts(), 4);
  
  std::vector<uint> n_prts_by_patch(mprts.n_patches());
  mprts.get_size_all(n_prts_by_patch.data());
  EXPECT_EQ(n_prts_by_patch, std::vector<uint>({n_prts, n_prts, n_prts, n_prts}));

  // check internal representation
  nn = 0;
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto& patch = mprts.grid().patches[p];
    for (auto prt: mprts.get_particles(p)) {
      // xm is patch-relative position
      auto xm = .5 * (patch.xe - patch.xb);
      EXPECT_EQ(prt.x[0], xm[0]);
      EXPECT_EQ(prt.x[1], xm[1]);
      EXPECT_EQ(prt.x[2], xm[2]);
      EXPECT_EQ(prt.w, nn);
      nn++;
    }
  }

  auto prts = mprts[0];
  auto cprts = mprts.get_particles(0, 1);
  typename Mparticles::patch_t::const_accessor prt{cprts[0], prts};

  auto prt2 = mprts[0].get_particle(0);
  EXPECT_EQ(cprts[0], prt2);
  
  nn = 0;
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto prts = mprts[p];
    auto& patch = mprts.grid().patches[p];

    for (auto prt: prts.get()) {
      auto x = .5 * (patch.xb + patch.xe);
      EXPECT_EQ(prt.position()[0], x[0]);
      EXPECT_EQ(prt.position()[1], x[1]);
      EXPECT_EQ(prt.position()[2], x[2]);
      EXPECT_EQ(prt.w(), nn);
      EXPECT_EQ(prt.kind(), 0);
      nn++;
    }
  }
}

#endif

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
