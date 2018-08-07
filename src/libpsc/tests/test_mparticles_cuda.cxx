
#include <gtest/gtest.h>

#include "psc_config.h"
#include "test_common.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/psc_particles_cuda.h"
#endif

template<typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using MparticlesCudaTestTypes = ::testing::Types<Config<MparticlesCuda<BS144>>>;

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
  const int n_prts = 1;
  
  auto mprts = this->mk_mprts();

  int nn = 0;
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto prts = mprts[p];
    auto& patch = mprts.grid().patches[p];
    for (int n = 0; n < n_prts; n++) {
      particle_inject prt = {};
      auto x = .5 * (patch.xb + patch.xe);
      int kind = 0;
      // use weight to store particle number for testing
      mprts.inject(p, {{x[0], x[1], x[2]}, {}, double(nn), kind});
      nn++;
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
