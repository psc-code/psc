
#include "cuda_mparticles.h"
#include "cuda_test.hxx"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <mrc_profile.h>

#include "gtest/gtest.h"

struct prof_globals prof_globals; // FIXME

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 0;
}

// ----------------------------------------------------------------------
// cuda_mparticles_add_particles_test_1
//
// add 1 particle at the center of each cell, in the "wrong" order in each
// patch (C order, but to get them ordered by block, they need to be reordered
// into Fortran order, a.k.a., this will exercise the initial sorting

struct SetParticleTest1
{
  SetParticleTest1(const Grid_t& grid)
  : grid_(grid)
  {
  }
  
  cuda_mparticles_prt operator()(int n)
  {
    Int3 ldims = grid_.ldims;
    Vec3<double> dx = grid_.dx;
    
    int k = n % ldims[2];
    n /= ldims[2];
    int j = n % ldims[1];
    n /= ldims[1];
    int i = n;

    cuda_mparticles_prt prt;
    prt.xi[0] = dx[0] * (i + .5f);
    prt.xi[1] = dx[1] * (j + .5f);
    prt.xi[2] = dx[2] * (k + .5f);
    prt.pxi[0] = i;
    prt.pxi[1] = j;
    prt.pxi[2] = k;
    prt.kind = 0;
    prt.qni_wni = 1.;
    return prt;
  }

private:
  const Grid_t& grid_;
};

void cuda_mparticles_add_particles_test_1(cuda_mparticles* cmprts, uint *n_prts_by_patch)
{
  const Grid_t& grid = cmprts->grid_;
  Int3 ldims = grid.ldims;

  for (int p = 0; p < cmprts->n_patches; p++) {
    n_prts_by_patch[p] = ldims[0] * ldims[1] * ldims[2];
  }

  SetParticleTest1 set_particles(grid);
  
  cmprts->reserve_all(n_prts_by_patch);
  cmprts->resize_all(n_prts_by_patch);
  
  for (int p = 0; p < grid.patches.size(); p++) {
    cmprts->set_particles(p, set_particles);
  }
}

// ----------------------------------------------------------------------
// get_particles_test

struct GetParticlesTest1
{
  void operator()(int n, cuda_mparticles_prt& prt) {
    printf("prt[%d] xi %g %g %g // pxi %g %g %g // kind %d // qni_wni %g\n",
	   n, prt.xi[0], prt.xi[1], prt.xi[2],
	   prt.pxi[0], prt.pxi[1], prt.pxi[2],
	   prt.kind, prt.qni_wni);
  }
};

void
get_particles_test(cuda_mparticles* cmprts)
{
  GetParticlesTest1 get_particles;
  for (int p = 0; p < cmprts->n_patches; p++) {
    cmprts->get_particles(p, get_particles);
  }
}

// ======================================================================
// CudaMparticlesTest

struct CudaMparticlesTest : TestBase, ::testing::Test
{
  std::unique_ptr<Grid_t> grid_;

  const Int3 bs_ = { 1, 1, 1 };

  void SetUp()
  {
    grid_.reset(new Grid_t({ 1, 4, 2 }, { 1., 40., 20. }));
  }
};

// ----------------------------------------------------------------------
TEST_F(CudaMparticlesTest, ConstructorDestructor)
{
  grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
  grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));
  std::unique_ptr<cuda_mparticles> cmprts(make_cmprts(*grid_));
  EXPECT_EQ(cmprts->n_patches, 1);
}

// ----------------------------------------------------------------------
TEST_F(CudaMparticlesTest, SetParticles)
{
  grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
  grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));
  std::unique_ptr<cuda_mparticles> cmprts(make_cmprts(*grid_));

  uint n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_add_particles_test_1(cmprts.get(), n_prts_by_patch);

  // check that particles are in C order
  cmprts->get_particles(0, [&] (int n, const cuda_mparticles_prt &prt) {
      int k = n % grid_->ldims[2]; n /= grid_->ldims[2];
      int j = n % grid_->ldims[1]; n /= grid_->ldims[1];
      int i = n;
      EXPECT_FLOAT_EQ(prt.xi[0], (i + .5) * grid_->dx[0]);
      EXPECT_FLOAT_EQ(prt.xi[1], (j + .5) * grid_->dx[1]);
      EXPECT_FLOAT_EQ(prt.xi[2], (k + .5) * grid_->dx[2]);
    });
}

// ----------------------------------------------------------------------
TEST_F(CudaMparticlesTest, SetupInternals)
{
  grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
  grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));
  std::unique_ptr<cuda_mparticles> cmprts(make_cmprts(*grid_));
  EXPECT_EQ(cmprts->n_patches, 1);

  uint n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_add_particles_test_1(cmprts.get(), n_prts_by_patch);

  // cmprts->dump();
  cmprts->setup_internals();
  // cmprts->dump();
  
  // check that particles are now in Fortran order
  cmprts->get_particles(0, [&] (int n, const cuda_mparticles_prt &prt) {
      int i = n % grid_->ldims[0]; n /= grid_->ldims[0];
      int j = n % grid_->ldims[1]; n /= grid_->ldims[1];
      int k = n;
      EXPECT_FLOAT_EQ(prt.xi[0], (i + .5) * grid_->dx[0]);
      EXPECT_FLOAT_EQ(prt.xi[1], (j + .5) * grid_->dx[1]);
      EXPECT_FLOAT_EQ(prt.xi[2], (k + .5) * grid_->dx[2]);
    });

  // FIXME, we should also check that h_off has been set correctly,
  // but really all of this should just call a consistency check
  // within cuda_mparticles
  // cmprts->dump();
}

// ---------------------------------------------------------------------
// SetupInternalsPieces
//
// Tests the pieces that go into setup_internals()

TEST_F(CudaMparticlesTest, SetupInternalsPieces)
{
  grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
  grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));

  std::vector<cuda_mparticles_prt> prts = {
    { .5, 35., 15. },
    { .5,  5.,  5. },
  };
  uint n_prts_by_patch[1];
  n_prts_by_patch[0] = prts.size();
  
  auto cmprts = make_cmprts(*grid_);
  cmprts->reserve_all(n_prts_by_patch);
  cmprts->resize_all(n_prts_by_patch);
  cmprts->set_particles(0, [&](int n) {
      return prts[n];
    });

  EXPECT_EQ(cmprts->d_bidx[0], 0);
  EXPECT_EQ(cmprts->d_bidx[1], 0);
  EXPECT_EQ(cmprts->d_id[0], 0);
  EXPECT_EQ(cmprts->d_id[1], 0);
  
  cmprts->find_block_indices_ids();
  
  EXPECT_EQ(cmprts->d_bidx[0], 7);
  EXPECT_EQ(cmprts->d_bidx[1], 0);
  EXPECT_EQ(cmprts->d_id[0], 0);
  EXPECT_EQ(cmprts->d_id[1], 1);

  cmprts->stable_sort_by_key();
  
  EXPECT_EQ(cmprts->d_bidx[0], 0);
  EXPECT_EQ(cmprts->d_bidx[1], 7);
  EXPECT_EQ(cmprts->d_id[0], 1);
  EXPECT_EQ(cmprts->d_id[1], 0);

  cmprts->reorder_and_offsets();

  float4 xi4_0 = cmprts->d_xi4[0], xi4_1 = cmprts->d_xi4[1];
  EXPECT_FLOAT_EQ(xi4_0.y, 5.);
  EXPECT_FLOAT_EQ(xi4_0.z, 5.);
  EXPECT_FLOAT_EQ(xi4_1.y, 35.);
  EXPECT_FLOAT_EQ(xi4_1.z, 15.);

  EXPECT_EQ(cmprts->d_off[0], 0);
  for (int b = 1; b < 8; b++) {
    EXPECT_EQ(cmprts->d_off[b], 1) << "where b = " << b;
  }
  EXPECT_EQ(cmprts->d_off[8], 2);
}

