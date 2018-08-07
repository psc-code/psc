
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

  void inject(Mparticles& mprts, int n_prts)
  {
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto prts = mprts[p];
      auto& patch = mprts.grid().patches[p];
      for (int n = 0; n < n_prts; n++) {
	double nn = double(n) / n_prts;
	auto L = patch.xe - patch.xb;
	particle_inject prt = {};
	prt.x[0] = patch.xb[0] + nn * L[0];
	prt.x[1] = patch.xb[1] + nn * L[1];
	prt.x[2] = patch.xb[2] + nn * L[2];
	prt.kind = 0;
	prt.w = 1.;
	mprts.inject(p, prt);
      }
    }
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


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
