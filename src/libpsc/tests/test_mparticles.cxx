
#include <gtest/gtest.h>

#define VPIC

#include "test_common.hxx"

#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_particles_vpic.h"

template<typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using MparticlesTestTypes = ::testing::Types<Config<MparticlesSingle>
					    ,Config<MparticlesDouble>
#ifdef VPIC
					    ,Config<MparticlesVpic>
#endif
					     >;

TYPED_TEST_CASE(MparticlesTest, MparticlesTestTypes);

// ======================================================================
// MparticlesTest

template<typename T>
struct MparticlesTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using MakeGrid = typename T::MakeGrid;
  
  Grid_t mk_grid()
  {
    Grid_t grid = MakeGrid{}();
    grid.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));

    return grid;
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
  
};

// -----------------------------------------------------------------------
// Constructor

TYPED_TEST(MparticlesTest, Constructor)
{
  using Mparticles = typename TypeParam::Mparticles;

  auto grid = this->mk_grid();
  Mparticles mprts(grid);
}

// ----------------------------------------------------------------------
// inject

TYPED_TEST(MparticlesTest, inject)
{
  using Mparticles = typename TypeParam::Mparticles;
  const int n_prts = 4;

  auto grid = this->mk_grid();
  Mparticles mprts(grid);

  this->inject(mprts, n_prts);
}

#ifndef VPIC // FIXME
// ----------------------------------------------------------------------
// setParticles

TYPED_TEST(MparticlesTest, setParticles)
{
  using Mparticles = typename TypeParam::Mparticles;
  const int n_prts = 4;

  auto grid = this->mk_grid();
  Mparticles mprts(grid);

  this->inject(mprts, n_prts);

  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto prts = mprts[p];
    auto& patch = mprts.grid().patches[p];
    EXPECT_EQ(prts.size(), n_prts);
    for (int n = 0; n < n_prts; n++) {
      double nn = double(n) / n_prts;
      auto L = patch.xe - patch.xb;
      auto prt = prts.get(n);
      auto x = prt.position();
      EXPECT_EQ(x[0], patch.xb[0] + nn * L[0]);
      EXPECT_EQ(x[1], patch.xb[1] + nn * L[1]);
      EXPECT_EQ(x[2], patch.xb[2] + nn * L[2]);
      EXPECT_EQ(prt.w(), 1.);
      EXPECT_EQ(prt.kind(), 0);
    }
  }
}
#endif

