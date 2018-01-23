
#include "cuda_mparticles.h"

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

struct CudaMparticlesTest : ::testing::Test
{
  std::unique_ptr<Grid_t> grid_;

  const Int3 bs_ = { 1, 1, 1 };

  void SetUp()
  {
    grid_.reset(new Grid_t({ 1, 4, 2 }, { 1., 40., 20. }));
  }

  cuda_mparticles* make_cmprts()
  {
    grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
    grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));
    struct cuda_mparticles *cmprts = new cuda_mparticles(*grid_, bs_);
    
    return cmprts;
  }
};

// ----------------------------------------------------------------------
TEST_F(CudaMparticlesTest, ConstructorDestructor)
{
  std::unique_ptr<cuda_mparticles> cmprts(make_cmprts());
  EXPECT_EQ(cmprts->n_patches, 1);
}

// ----------------------------------------------------------------------
TEST_F(CudaMparticlesTest, FullInitialSetup)
{
  std::unique_ptr<cuda_mparticles> cmprts(make_cmprts());

  uint n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_add_particles_test_1(cmprts.get(), n_prts_by_patch);
  printf("added particles\n");
  cmprts->dump_by_patch(n_prts_by_patch);

  cmprts->setup_internals();
  printf("set up internals\n");
  cmprts->dump();

  printf("get_particles_test\n");
  get_particles_test(cmprts.get());
}
