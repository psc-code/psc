
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

TEST_F(CudaMparticlesBndTest, BndPrep)
{
  grid_->kinds.push_back(Grid_t::Kind(1., 1., "test species"));

  std::unique_ptr<cuda_mparticles> cmprts(make_cmprts(*grid_));

  // (ab)use kind to track particle more easily in the test
  std::vector<cuda_mparticles_prt> prts = {
    { .5,  5., 5., 0., 0., 0., 0 },
    { .5, 35., 5., 0., 0., 0., 1 },

    { .5,  5., 5., 0., 0., 0., 2 },
    { .5, 35., 5., 0., 0., 0., 3 },
  };

  uint n_prts_by_patch[cmprts->n_patches];
  n_prts_by_patch[0] = 2;
  n_prts_by_patch[1] = 2;

  // FIXME eventually shouldn't have to reserve additional room for sending here
  uint n_prts_reserve_by_patch[cmprts->n_patches];
  n_prts_reserve_by_patch[0] = 2;
  n_prts_reserve_by_patch[1] = 4;
  
  cmprts->reserve_all(n_prts_reserve_by_patch);
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

  // test spine_reduce
  cmprts->spine_reduce(cmprts.get());

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

  // test find_n_send
  cmprts->find_n_send(cmprts.get());

  for (int p = 0; p < cmprts->n_patches; p++) {
    printf("p %d: n_send %d\n", p, cmprts->bpatch[p].n_send);
    EXPECT_EQ(cmprts->bpatch[p].n_send, 1);
  }
  EXPECT_EQ(cmprts->n_prts_send, 2);

  // test scan_send_buf_total
  cmprts->scan_send_buf_total_gold(cmprts.get());

#if 1
  printf("sums: ");
  for (int n = 0; n < cmprts->n_prts; n++) {
    int sum = cmprts->d_sums[n];
    printf(" %3d", sum);
  }
  printf("\n");
#endif

  // where in the send region at the tail the OOB particles should go
  EXPECT_EQ(cmprts->d_sums[1], 0);
  EXPECT_EQ(cmprts->d_sums[3], 1);

  // particles 1, 3, which need to be exchanged, should now be at the
  // end of the regular array
  EXPECT_EQ(cuda_float_as_int(float4(cmprts->d_xi4[cmprts->n_prts  ]).w), 1);
  EXPECT_EQ(cuda_float_as_int(float4(cmprts->d_xi4[cmprts->n_prts+1]).w), 3);

  // test copy_from_dev_and_convert
  cmprts->copy_from_dev_and_convert(cmprts.get());

#if 0
  for (int p = 0; p < cmprts->n_patches; p++) {
    printf("from_dev: p %d\n", p);
    for (auto& prt : cmprts->bpatch[p].buf) {
      printf("  prt xyz %g %g %g kind %d\n", prt.xi, prt.yi, prt.zi, prt.kind_);
    }
  }
#endif

  EXPECT_EQ(cmprts->bpatch[0].buf.size(), 1);
  EXPECT_EQ(cmprts->bpatch[1].buf.size(), 1);
  EXPECT_EQ(cmprts->bpatch[0].buf[0].kind_, 1);
  EXPECT_EQ(cmprts->bpatch[1].buf[0].kind_, 3);
}