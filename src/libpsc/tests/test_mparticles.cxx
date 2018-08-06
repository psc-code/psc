
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

  MparticlesTest()
    : grid_{MakeGrid{}()}
  {
    grid_.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));

    // Setup basic grid parameters
    auto& domain = grid_.domain;
    double dx[3], xl[3], xh[3];
    for (int d = 0; d < 3; d++) {
      dx[d] = domain.length[d] / domain.gdims[d];
      xl[d] = domain.corner[d];
      xh[d] = xl[d] + domain.length[d];
    }

    vgrid_.setup(dx, grid_.dt, 1., 1.);
  
    // Define the grid
    //sim->define_periodic_grid(xl, xh, domain.gdims, domain.np);
  }

  template<typename tag>
  Mparticles mk_mprts(tag dummy)
  {
    Mparticles mprts(grid_);
    mprts.define_species("test_species", 1., 1., 100, 10,
			 10, 0);
    return mprts;
  }

  Mparticles mk_mprts(MparticlesVpic* dummy)
  {
    Mparticles mprts(grid_, &vgrid_);
    mprts.define_species("test_species", 1., 1., 100, 10,
			 10, 0);
    return mprts;
  }

  Mparticles mk_mprts()
  {
    return mk_mprts(static_cast<Mparticles*>(nullptr));
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
  Grid vgrid_;
};

// -----------------------------------------------------------------------
// Constructor

TYPED_TEST(MparticlesTest, Constructor)
{
  auto mprts = this->mk_mprts();
}

// ----------------------------------------------------------------------
// inject

TYPED_TEST(MparticlesTest, inject)
{
  const int n_prts = 4;

  auto mprts = this->mk_mprts();

  this->inject(mprts, n_prts);
}

#ifndef VPIC // FIXME
// ----------------------------------------------------------------------
// setParticles

TYPED_TEST(MparticlesTest, setParticles)
{
  const int n_prts = 4;

  auto mprts = this->mk_mprts();

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

