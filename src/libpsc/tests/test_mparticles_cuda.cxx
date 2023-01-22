
#include <gtest/gtest.h>

#include <mpi.h>

#ifdef USE_CUDA

#include "test_common.hxx"
#include "PscConfig.h"
#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "cuda_base.hxx"

// FIXME, the general tests should be moved -> test_mparticles,
// and the real cuda ones to test_cuda_mparticle?

#include "psc_particles_single.h"

template <typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using MparticlesCudaTestTypes =
  ::testing::Types<Config<MparticlesCuda<BS144>, MakeTestGridYZ>>;

TYPED_TEST_SUITE(MparticlesCudaTest, MparticlesCudaTestTypes);

// ======================================================================
// MparticlesCudaTest

template <typename T>
struct MparticlesCudaTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using MakeGrid = typename T::MakeGrid;

  MparticlesCudaTest() : grid_{MakeGrid{}()}
  {
    grid_.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  }

  Mparticles mk_mprts()
  {
    Mparticles mprts(grid_);
    mprts.define_species("test_species", 1., 1., 100, 10, 10, 0);
    return mprts;
  }

  const Grid_t& grid() const { return grid_; }

private:
  Grid_t grid_;
};

// ----------------------------------------------------------------------
// Conversion from MparticlesSingle

TYPED_TEST(MparticlesCudaTest, ConvertFromSingle)
{
  using Mparticles = typename TypeParam::Mparticles;
  const int n_prts = 1;
  const auto& grid = this->grid();

  auto mprts_single = MparticlesSingle{grid};

  EXPECT_EQ(mprts_single.n_patches(), 4);

  int nn = 0;
  {
    auto inj = mprts_single.injector();
    for (int p = 0; p < mprts_single.n_patches(); ++p) {
      auto injector = inj[p];
      auto& patch = grid.patches[p];
      for (int n = 0; n < n_prts; n++) {
        auto x = .5 * (patch.xb + patch.xe);
        int kind = 0;
        // use weight to store particle number for testing
        injector({{x[0], x[1], x[2]}, {}, double(nn), kind});
        nn++;
      }
    }
  }

  EXPECT_EQ(mprts_single.size(), 4);

  auto mprts = mprts_single.get_as<Mparticles>();

  EXPECT_EQ(mprts.size(), 4);

  nn = 0;
  auto accessor = mprts.accessor();
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto prts = accessor[p];
    auto& patch = mprts.grid().patches[p];

    for (auto prt : prts) {
      auto x = .5 * (patch.xb + patch.xe);
      EXPECT_EQ(prt.position()[0], x[0]);
      EXPECT_EQ(prt.position()[1], x[1]);
      EXPECT_EQ(prt.position()[2], x[2]);
      EXPECT_EQ(prt.w(), nn);
      EXPECT_EQ(prt.kind(), 0);
      nn++;
    }
  }

  auto&& mprts2 = mprts.template get_as<MparticlesSingle>();
  {
    nn = 0;
    auto accessor = mprts2.accessor();
    for (int p = 0; p < mprts2.n_patches(); ++p) {
      auto prts = accessor[p];
      auto& patch = mprts2.grid().patches[p];

      for (auto prt : prts) {
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
}

#endif

// ======================================================================
// main

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
#ifdef USE_CUDA
  cuda_base_init();
#endif
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
