
#include "cuda_mparticles.h"
#include "cuda_test.hxx"

#include <mrc_profile.h>

#include "gtest/gtest.h"

struct prof_globals prof_globals; // FIXME

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 0;
}

// ======================================================================
// CudaMparticlesBndTest

struct CudaMparticlesBndTest : TestBase, ::testing::Test
{
  using Double3 = Vec3<double>;
  
  std::unique_ptr<Grid_t> grid_;

  const Int3 bs_ = { 1, 1, 1 };

  void SetUp()
  {
    grid_.reset(new Grid_t());
    grid_->gdims = { 1, 8, 8 };
    grid_->ldims = { 1, 4, 4 };
    grid_->dx = Double3{ 1., 80., 80. } / Double3(grid_->gdims);
    grid_->patches.emplace_back(Grid_t::Patch({ 0.,  0.,  0.}, { 1., 40., 40. }));
    grid_->patches.emplace_back(Grid_t::Patch({ 0., 40.,  0.}, { 1., 80., 40. }));
  }
};

TEST_F(CudaMparticlesBndTest, SpineReduce)
{
  grid_->kinds.push_back(Grid_t::Kind(1., 1., "test species"));

  std::unique_ptr<cuda_mparticles> cmprts(make_cmprts(*grid_));

  std::vector<cuda_mparticles_prt> prts = {
    { .5,  5., 5. },
    { .5, 35., 5. },

    { .5,  5., 5. },
    { .5, 35., 5. },
  };

  uint n_prts_by_patch[cmprts->n_patches];
  n_prts_by_patch[0] = 2;
  n_prts_by_patch[1] = 2;
  
  cmprts->reserve_all(n_prts_by_patch);
  cmprts->inject(prts.data(), n_prts_by_patch);

  // move every particle one full cell to the right (+y, that is)
  // (position doesn't actually matter since we'll only look at bidx)
  for (int n = 0; n < cmprts->n_prts; n++) {
    float4 xi4 = cmprts->d_xi4[n];
    xi4.y += 10.;
    cmprts->d_xi4[n] = xi4;
  }
  cmprts->d_bidx[0] = 0 + 1 * 3; // +1 in y, 0 in z
  cmprts->d_bidx[1] = CUDA_BND_S_OOB;
  cmprts->d_bidx[2] = 0 + 1 * 3; // +1 in y, 0 in z
  cmprts->d_bidx[3] = CUDA_BND_S_OOB;

#if 0
  cmprts->dump();
#endif

  cmprts->spine_reduce_gold(cmprts.get());

#if 0
  for (int b = 0; b < cmprts->n_blocks; b++) {
    printf("b %d:", b);
    for (int n = 0; n < 10; n++) {
      int cnt = cmprts->d_spine_cnts[10*b + n];
      printf(" %3d", cnt);
    }
    printf("\n");
  }
#endif

  for (int b = 0; b < cmprts->n_blocks; b++) {
    for (int n = 0; n < 10; n++) {
      int cnt = cmprts->d_spine_cnts[10*b + n];
      // one particle each moves to block 1, 17, respectively, from the left (-y: 3)
      if ((b ==  1 && n == 3) ||
	  (b == 17 && n == 3)) {
	EXPECT_EQ(cnt, 1) << "where b = " << b << " n = " << n;
      } else {
	EXPECT_EQ(cnt, 0) << "where b = " << b << " n = " << n;
      }
    }
  }

#if 0
  printf("oob: ");
  for (int b = 0; b < cmprts->n_blocks + 1; b++) {
    int cnt = cmprts->d_spine_cnts[10*cmprts->n_blocks + b];
    printf(" %3d", cnt);
  }
  printf("\n");
#endif

  for (int b = 0; b < cmprts->n_blocks + 1; b++) {
    int cnt = cmprts->d_spine_cnts[10*cmprts->n_blocks + b];
    // the particles in cell 3 and 19 went out of bounds
    if (b == 3 || b == 19) {
      EXPECT_EQ(cnt, 1) << "where b = " << b;
    } else {
      EXPECT_EQ(cnt, 0) << "where b = " << b;
    }
  }

#if 0
  printf("sum: ");
  for (int b = 0; b < cmprts->n_blocks + 1; b++) {
    int cnt = cmprts->d_spine_sums[10*cmprts->n_blocks + b];
    printf(" %3d", cnt);
  }
  printf("\n");
#endif
  
  for (int b = 0; b < cmprts->n_blocks + 1; b++) {
    int cnt = cmprts->d_spine_sums[10*cmprts->n_blocks + b];
    // the particles in cell 3 and 19 went out of bounds
    if (b <= 3) {
      EXPECT_EQ(cnt, 0) << "where b = " << b;
    } else if (b <= 19) {
      EXPECT_EQ(cnt, 1) << "where b = " << b;
    } else {
      EXPECT_EQ(cnt, 2) << "where b = " << b;
    }
  }
}