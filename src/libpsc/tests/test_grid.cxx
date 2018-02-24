
#include <gtest/gtest.h>

#include "grid.hxx"

static Grid_t make_grid()
{
  Int3 gdims = { 8, 4, 2 };
  Int3 ldims = { 4, 2, 2 };
  Grid_t::Real3 length = {  80.,  40., 20. };
  Grid_t::Real3 corner = { -40., -20., 0. };
  std::vector<Int3> offs = { { 0, 0, 0 }, { 4, 0, 0 } };
  return Grid_t(gdims, ldims, length, corner, offs);
}

static Grid_t make_grid_1()
{
  Int3 gdims = { 1, 1, 1 };
  Int3 ldims = { 1, 1, 1 };
  Grid_t::Real3 length = {  10., 10., 10. };
  Grid_t::Real3 corner = {   0.,  0.,  0. };
  std::vector<Int3> offs = { { 0, 0, 0 } };
  return Grid_t(gdims, ldims, length, corner, offs);
}

TEST(Grid, Domain)
{
  Grid_t grid = make_grid();
    
  EXPECT_EQ(grid.gdims, Int3({ 8, 4, 2 }));
  EXPECT_EQ(grid.ldims, Int3({ 4, 2, 2 }));
  EXPECT_EQ(grid.dx, Grid_t::Real3({ 10., 10., 10. }));
  EXPECT_EQ(grid.n_patches(), 2);
  EXPECT_EQ(grid.patches[0].off, Int3({ 0, 0, 0 }));
  EXPECT_EQ(grid.patches[0].xb, Grid_t::Real3({ -40., -20.,  0. }));
  EXPECT_EQ(grid.patches[0].xe, Grid_t::Real3({   0.,   0., 20. }));
  EXPECT_EQ(grid.patches[1].off, Int3({ 4, 0, 0 }));
  EXPECT_EQ(grid.patches[1].xb, Grid_t::Real3({   0., -20.,  0. }));
  EXPECT_EQ(grid.patches[1].xe, Grid_t::Real3({  40.,   0., 20. }));
}

TEST(Grid, Kinds)
{
  Grid_t grid = make_grid();

  grid.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  EXPECT_EQ(grid.kinds.size(), 1);
}

#include "psc.h" // FIXME, just for EX, ...
#include "psc_fields_single.h"

template<typename mfields, typename Set>
void setValues(mfields& mflds, Set set)
{
  for (int p = 0; p < mflds.n_patches(); p++) {
    auto F = mflds[p];
    F(EX, 0,0,0) = set(EX);
    F(EX, 0,1,0) = set(EX);
    F(EX, 0,0,1) = set(EX);
    F(EX, 0,1,1) = set(EX);

    F(EY, 0,0,0) = set(EY);
    F(EY, 0,0,1) = set(EY);
    F(EY, 1,0,0) = set(EY);
    F(EY, 1,0,1) = set(EY);
    
    F(EZ, 0,0,0) = set(EZ);
    F(EZ, 1,0,0) = set(EZ);
    F(EZ, 0,1,0) = set(EZ);
    F(EZ, 1,1,0) = set(EZ);

    F(HX, 0,0,0) = set(HX);
    F(HX, 1,0,0) = set(HX);

    F(HY, 0,0,0) = set(HY);
    F(HY, 0,1,0) = set(HY);

    F(HZ, 0,0,0) = set(HZ);
    F(HZ, 0,0,1) = set(HZ);
  }
}
  
TEST(mfields, Constructor)
{
  using Mfields = psc_mfields_<fields_single_t>;

  Grid_t grid = make_grid();
  Mfields mflds(grid, NR_FIELDS, Int3{ 1, 1, 1 });

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());
}

TEST(mfields, Set)
{
  using Mfields = psc_mfields_<fields_single_t>;

  Grid_t grid = make_grid();
  Mfields mflds(grid, NR_FIELDS, Int3{ 1, 1, 1 });

  setValues(mflds, [](int m) -> Mfields::real_t {
      switch(m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  auto F = mflds[0];
  EXPECT_EQ(F(EY, 0, 0, 0), 2.);
}

#include "../vpic/PscRng.h"

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

TEST(Rng, RngPool)
{
  RngPool rngpool;
  Rng *rng = rngpool[0];

  for (int i = 0; i < 100; ++i) {
    double r = rng->uniform(0., 1000.);
    EXPECT_GE(r, 0.);
    EXPECT_LE(r, 1000.);
  }
}

#include "psc_particles_single.h"

template<typename Mparticles>
static void initMparticlesRandom(Mparticles& mprts, int n_prts)
{
  RngPool rngpool;
  Rng *rng = rngpool[0];

  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto patch = mprts.grid().patches[p];
    for (int n = 0; n < n_prts; n++) {
      psc_particle_inject prt = {};
      prt.x[0] = rng->uniform(patch.xb[0], patch.xe[0]);
      prt.x[1] = rng->uniform(patch.xb[1], patch.xe[1]);
      prt.x[2] = rng->uniform(patch.xb[2], patch.xe[2]);
      prt.kind = 0;
      prt.w = 1.;
      mprts.inject(p, prt);
    }
  }
}

TEST(Mparticles, Constructor)
{
  using Mparticles = psc_mparticles_<particle_single_t>;
  
  Grid_t grid = make_grid();
  grid.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));

  Mparticles mprts(grid);
}

TEST(Mparticles, setParticles)
{
  using particle_t = particle_single_t;
  using Mparticles = psc_mparticles_<particle_t>;
  const int n_prts = 131;

  Grid_t grid = make_grid();
  grid.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));

  Mparticles mprts(grid);
  initMparticlesRandom(mprts, n_prts);

  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto& prts = mprts[p];
    EXPECT_EQ(prts.size(), n_prts);
    for (int n = 0; n < n_prts; n++) {
      particle_t& prt = prts[n];
      EXPECT_EQ(mprts.prt_wni(prt), 1.);
      EXPECT_EQ(prt.kind(), 0);
    }
  }
}

#include "../libpsc/psc_push_particles/1vb/push_particles_1vbec_single.hxx"

TEST(PushParticles, Accel)
{
  using Mparticles = psc_mparticles_<particle_single_t>;
  using Mfields = psc_mfields_<fields_single_t>;
  using PushParticles = PushParticles1vbecSingle;
  const int n_prts = 131;
  const int n_steps = 10;
  const Mparticles::real_t eps = 1e-6;
  
  Grid_t grid = make_grid_1();
  grid.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));

  Mfields mflds(grid, NR_FIELDS, Int3{ 1, 1, 1 });

  setValues(mflds, [](int m) -> Mfields::real_t {
      switch(m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  Mparticles mprts(grid);
  initMparticlesRandom(mprts, n_prts);
  
  PushParticles pushp;

  MPI_Init(0, 0);

  psc psc{};
  psc.grid_ = &grid;
  ppsc = &psc;
  
  psc_mparticles _mprts[1];
  _mprts->obj.subctx = &mprts;

  extern struct psc_mfields_ops psc_mfields_single_ops;
  psc_mfields _mflds[1];
  _mflds->obj.subctx = &mflds;
  _mflds->obj.ops = (mrc_obj_ops*) &psc_mfields_single_ops;
  
  // run test
  for (int n = 0; n < n_steps; n++) {
    pushp.push_mprts(_mprts, _mflds);
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto& prts = mprts[p];
      for (int m = 0; m < prts.size(); m++) {
	auto& prt = prts[m];
    	EXPECT_NEAR(prt.pxi, 1*(n+1), eps);
    	EXPECT_NEAR(prt.pyi, 2*(n+1), eps);
    	EXPECT_NEAR(prt.pzi, 3*(n+1), eps);
      };
    }
  }
}
