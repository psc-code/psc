
#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "../libpsc/vpic/mparticles_vpic.hxx"
#include "../libpsc/vpic/PscGridBase.h"
#include "../libpsc/vpic/PscParticleBc.h"
#include "../libpsc/vpic/PscParticlesBase.h"
#include "../libpsc/vpic/vpic_config.h"

template<typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using BalanceTestTypes = ::testing::Types<Config<MparticlesSingle>
					 ,Config<MparticlesSingle, MakeTestGridYZ>
					 ,Config<MparticlesDouble>
					 >;

TYPED_TEST_CASE(BalanceTest, BalanceTestTypes);

// ======================================================================
// BalanceTest

template<typename T>
struct BalanceTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using Particle = typename Mparticles::Particle;
  using MakeGrid = typename T::MakeGrid;

  BalanceTest()
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

  template<typename _Mparticles>
  void inject_test_particles(_Mparticles& mprts, int n_prts)
  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];
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
	injector(prt);
      }
    }
  }

  const Grid_t& grid() { return grid_; }
  
private:
  Grid_t grid_;
};

// -----------------------------------------------------------------------
// Constructor

TYPED_TEST(BalanceTest, Constructor)
{
  auto mprts = this->mk_mprts();
}

// ----------------------------------------------------------------------
// Inject

TYPED_TEST(BalanceTest, Inject)
{
  const int n_prts = 4;

  auto mprts = this->mk_mprts();

  this->inject_test_particles(mprts, n_prts);
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
