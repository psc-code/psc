
#include "cuda_mparticles.h"
#include "cuda_bndp.h"
#include "cuda_test.hxx"

#include <mrc_profile.h>

#include "gtest/gtest.h"

struct prof_globals prof_globals; // FIXME

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 0;
}

using CudaMparticles = cuda_mparticles<BS144>;

// ======================================================================
// CudaMparticlesBndTest

struct CudaMparticlesBndTest : TestBase<CudaMparticles>, ::testing::Test
{
  using Double3 = Vec3<double>;
  
  std::unique_ptr<Grid_t> grid;
  std::unique_ptr<CudaMparticles> cmprts;
  std::unique_ptr<cuda_bndp<CudaMparticles, dim_yz>> cbndp;

  void SetUp()
  {
    auto domain = Grid_t::Domain{{1, 32, 32}, {1., 320., 320.}, {0., 0., 0.},
				 {1, 2, 2}};
    auto bc = GridBc{};
    auto kinds = Grid_t::Kinds{Grid_t::Kind{1., 1., "k0"},
			       Grid_t::Kind{1., 1., "k1"},
			       Grid_t::Kind{1., 1., "k2"},
			       Grid_t::Kind{1., 1., "k3"}};
    auto norm = Grid_t::Normalization{};
    double dt = .1;
    grid.reset(new Grid_t(domain, bc, kinds, norm, dt));

    grid->kinds.push_back(Grid_t::Kind(1., 1., "test species"));

    cmprts.reset(make_cmprts(*grid));

    // (ab)use kind to track particle more easily in the test
    using Real3 = cuda_mparticles_prt::Real3;
    std::vector<cuda_mparticles_prt> prts = {
      {{ .5,  35., 5.}, {}, 0., 0},
      {{ .5, 155., 5.}, {}, 0., 1},
      
      {{ .5,  35., 5.}, {}, 0., 2},
      {{ .5, 155., 5.}, {}, 0., 3},
    };

    std::vector<uint> n_prts_by_patch = {2, 2, 0, 0};
    
    // FIXME eventually shouldn't have to reserve additional room for sending here
    std::vector<uint> n_prts_reserve_by_patch = {2, 4, 0, 0};
    
    cmprts->reserve_all(n_prts_reserve_by_patch.data());
    cmprts->inject_buf(prts.data(), n_prts_by_patch.data());

    // move every particle one full cell to the right (+y, that is)
    // (position doesn't actually matter since we'll only look at bidx)
    for (int n = 0; n < cmprts->n_prts; n++) {
      float4 xi4 = cmprts->d_xi4[n];
      xi4.y += 10.;
      cmprts->d_xi4[n] = xi4;
    }
    auto& d_bidx = cmprts->by_block_.d_idx;
    d_bidx[0] = 0 + 1 * 3; // +1 in y, 0 in z
    d_bidx[1] = CUDA_BND_S_OOB;
    d_bidx[2] = 0 + 1 * 3; // +1 in y, 0 in z
    d_bidx[3] = CUDA_BND_S_OOB;
    
#if 0
    cmprts->dump();
#endif

    cbndp.reset(new cuda_bndp<cuda_mparticles<BS144>, dim_yz>(*grid));
  }
};


// ----------------------------------------------------------------------
// BndPrep
//
// tests cuda_bndp::prep()

TEST_F(CudaMparticlesBndTest, BndPrep)
{
  cbndp->prep(cmprts.get());

  // particles 0 and 2 remain in their patch,
  // particles 1 and 3 leave their patch and need special handling
  EXPECT_EQ(cbndp->bpatch[0].buf.size(), 1);
  EXPECT_EQ(cbndp->bpatch[1].buf.size(), 1);
  EXPECT_EQ(cbndp->bpatch[0].buf[0].kind, 1);
  EXPECT_EQ(cbndp->bpatch[1].buf[0].kind, 3);
}

// ----------------------------------------------------------------------
// BndPrepDetail
//
// tests the pieces that go into cuda_bndp::prep()

TEST_F(CudaMparticlesBndTest, BndPrepDetail)
{
  // test spine_reduce
  cbndp->spine_reduce(cmprts.get());

#if 0
  for (int b = 0; b < cmprts->n_blocks; b++) {
    printf("b %d:", b);
    for (int n = 0; n < 10; n++) {
      int cnt = cbndp->d_spine_cnts[10*b + n];
      printf(" %3d", cnt);
    }
    printf("\n");
  }
#endif

  for (int b = 0; b < cmprts->n_blocks; b++) {
    for (int n = 0; n < 10; n++) {
      int cnt = cbndp->d_spine_cnts[10*b + n];
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
    int cnt = cbndp->d_spine_cnts[10*cmprts->n_blocks + b];
    printf(" %3d", cnt);
  }
  printf("\n");
#endif

  for (int b = 0; b < cmprts->n_blocks + 1; b++) {
    int cnt = cbndp->d_spine_cnts[10*cmprts->n_blocks + b];
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
    int cnt = cbndp->d_spine_sums[10*cmprts->n_blocks + b];
    printf(" %3d", cnt);
  }
  printf("\n");
#endif
  
  for (int b = 0; b < cmprts->n_blocks + 1; b++) {
    int cnt = cbndp->d_spine_sums[10*cmprts->n_blocks + b];
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
  cbndp->n_prts_send = cbndp->find_n_send(cmprts.get());

  for (int p = 0; p < cmprts->n_patches; p++) {
    //printf("p %d: n_send %d\n", p, cmprts->bpatch[p].n_send);
    EXPECT_EQ(cbndp->bpatch[p].n_send, p < 2 ? 1 : 0);
  }
  EXPECT_EQ(cbndp->n_prts_send, 2);

  // test scan_send_buf_total
#if 1
  cbndp->scan_send_buf_total(cmprts.get(), cbndp->n_prts_send);

#if 0
  printf("ids: ");
  for (int n = cmprts->n_prts - cmprts->n_prts_send; n < cmprts->n_prts; n++) {
    int id = cmprts->d_id[n];
    printf(" %3d", id);
  }
  printf("\n");
#endif
  EXPECT_EQ(cmprts->n_prts, 4);
  EXPECT_EQ(cbndp->n_prts_send, 2);
  EXPECT_EQ(cmprts->by_block_.d_id[2], 1);
  EXPECT_EQ(cmprts->by_block_.d_id[3], 3);

#else
  cbndp->scan_send_buf_total_gold(cmprts.get(), cbndp->n_prts_send);
  // the intermediate scan_send_buf_total_gold result
  // can be tested here, but the non-gold version works differently
  // and has different intermediate results
#if 0
  printf("sums: ");
  for (int n = 0; n < cmprts->n_prts; n++) {
    int sum = cmprts->d_sums[n];
    printf(" %3d", sum);
  }
  printf("\n");
#endif

  // where in the send region at the tail the OOB particles should go
  EXPECT_EQ(cbndp->d_sums[1], 0);
  EXPECT_EQ(cbndp->d_sums[3], 1);
#endif

  // particles 1, 3, which need to be exchanged, should now be at the
  // end of the regular array
  EXPECT_EQ(cuda_float_as_int(float4(cmprts->d_xi4[cmprts->n_prts  ]).w), 1);
  EXPECT_EQ(cuda_float_as_int(float4(cmprts->d_xi4[cmprts->n_prts+1]).w), 3);

  // test copy_from_dev_and_convert
  cbndp->copy_from_dev_and_convert(cmprts.get(), cbndp->n_prts_send);

#if 0
  for (int p = 0; p < cmprts->n_patches; p++) {
    printf("from_dev: p %d\n", p);
    for (auto& prt : cmprts->bpatch[p].buf) {
      printf("  prt xyz %g %g %g kind %d\n", prt.xi, prt.yi, prt.zi, prt.kind_);
    }
  }
#endif

  EXPECT_EQ(cbndp->bpatch[0].buf.size(), 1);
  EXPECT_EQ(cbndp->bpatch[1].buf.size(), 1);
  EXPECT_EQ(cbndp->bpatch[0].buf[0].kind, 1);
  EXPECT_EQ(cbndp->bpatch[1].buf[0].kind, 3);
}

// ----------------------------------------------------------------------
// BndPost
//
// tests cuda_bndp::post()

TEST_F(CudaMparticlesBndTest, BndPost)
{
  // BndPost expects the work done by bnd_prep()
  cbndp->prep(cmprts.get());

  // particles 0 and 2 remain in their patch,
  // particles 1 and 3 leave their patch and need special handling
  EXPECT_EQ(cbndp->bpatch[0].buf.size(), 1);
  EXPECT_EQ(cbndp->bpatch[1].buf.size(), 1);
  EXPECT_EQ(cbndp->bpatch[0].buf[0].kind, 1);
  EXPECT_EQ(cbndp->bpatch[1].buf[0].kind, 3);

  // Mock what the actual boundary exchange does, ie., move
  // particles to their new patch and adjust the relative position.
  // This assumes periodic b.c.
  particle_cuda_t prt1 = cbndp->bpatch[0].buf[0];
  particle_cuda_t prt3 = cbndp->bpatch[1].buf[0];
  prt1.x[1] -= 40.;
  prt3.x[1] -= 40.;
  cbndp->bpatch[0].buf[0] = prt3;
  cbndp->bpatch[1].buf[0] = prt1;
  
  cbndp->post(cmprts.get());

  // bnd_post doesn't do the actual final reordering
  EXPECT_TRUE(cmprts->need_reorder);
  cmprts->reorder();
  EXPECT_TRUE(cmprts->check_ordered());

#if 0
  cmprts->dump();
#endif
}

// ----------------------------------------------------------------------
// BndPostDetail
//
// tests the pieces that go into cuda_bndp::post()

TEST_F(CudaMparticlesBndTest, BndPostDetail)
{
  // BndPost expects the work done by bnd_prep()
  cbndp->prep(cmprts.get());

  // particles 0 and 2 remain in their patch,
  // particles 1 and 3 leave their patch and need special handling
  EXPECT_EQ(cbndp->bpatch[0].buf.size(), 1);
  EXPECT_EQ(cbndp->bpatch[1].buf.size(), 1);
  EXPECT_EQ(cbndp->bpatch[0].buf[0].kind, 1);
  EXPECT_EQ(cbndp->bpatch[1].buf[0].kind, 3);

  // Mock what the actual boundary exchange does, ie., move
  // particles to their new patch and adjust the relative position.
  // This assumes periodic b.c.
  particle_cuda_t prt1 = cbndp->bpatch[0].buf[0];
  particle_cuda_t prt3 = cbndp->bpatch[1].buf[0];
  prt1.x[1] -= 160.;
  prt3.x[1] -= 160.;
  cbndp->bpatch[0].buf[0] = prt3;
  cbndp->bpatch[1].buf[0] = prt1;

  // === test convert_and_copy_to_dev()
  uint n_prts_recv = cbndp->convert_and_copy_to_dev(cmprts.get());
  cmprts->n_prts += n_prts_recv;

  // n_recv should be set for each patch, and its total
  EXPECT_EQ(cbndp->bpatch[0].n_recv, 1);
  EXPECT_EQ(cbndp->bpatch[1].n_recv, 1);
  EXPECT_EQ(n_prts_recv, 2);

  // the received particle have been added to the previous total
  EXPECT_EQ(cmprts->n_prts, 6);

  // and the particle have been appended after the old end of the particle list
  int n_prts_old = cmprts->n_prts - n_prts_recv;
  EXPECT_EQ(cuda_float_as_int(float4(cmprts->d_xi4[n_prts_old  ]).w), 3);
  EXPECT_EQ(cuda_float_as_int(float4(cmprts->d_xi4[n_prts_old+1]).w), 1);

  // block indices have been calculated
  auto& d_bidx = cmprts->by_block_.d_idx;
  EXPECT_EQ(d_bidx[n_prts_old  ], 0);  // 0th block in 0th patch
  EXPECT_EQ(d_bidx[n_prts_old+1], 16); // 0th block in 1st patch

  // received particles per block have been counted
  for (int b = 0; b < cmprts->n_blocks; b++) {
    if (b == 0 || b == 16) {
      EXPECT_EQ(cbndp->d_spine_cnts[10*cmprts->n_blocks + b], 1);
    } else {
      EXPECT_EQ(cbndp->d_spine_cnts[10*cmprts->n_blocks + b], 0);
    }
  }

  // both particles are the 0th (and only) particle added to their respective block
  EXPECT_EQ(cbndp->d_bnd_off[0], 0);
  EXPECT_EQ(cbndp->d_bnd_off[1], 0);

  // === test sort
  uint n_prts_by_patch[cmprts->n_patches];
  cmprts->get_size_all(n_prts_by_patch);
  EXPECT_EQ(n_prts_by_patch[0], 2);
  EXPECT_EQ(n_prts_by_patch[1], 2);
  
  cbndp->sort_pairs_device(cmprts.get(), n_prts_recv);
  cmprts->n_prts -= cbndp->n_prts_send;

  EXPECT_EQ(cmprts->n_prts, 4);
  auto& d_id = cmprts->by_block_.d_id;
  EXPECT_EQ(d_id[0], 4);
  EXPECT_EQ(d_id[1], 0);
  EXPECT_EQ(d_id[2], 5);
  EXPECT_EQ(d_id[3], 2);

  cbndp->update_offsets(cmprts.get());
  auto& d_off = cmprts->by_block_.d_off;
  for (int b = 0; b <= cmprts->n_blocks; b++) {
    //if (b < cmprts->n_blocks) printf("b %d: off [%d:%d[\n", b, int(d_off[b]), int(d_off[b+1]));
    if (b < 1) {
      EXPECT_EQ(d_off[b], 0) << "where b = " << b;
    } else if (b < 2) {
      EXPECT_EQ(d_off[b], 1) << "where b = " << b;
    } else if (b < 17) {
      EXPECT_EQ(d_off[b], 2) << "where b = " << b;
    } else if (b < 18) {
      EXPECT_EQ(d_off[b], 3) << "where b = " << b;
    } else {
      EXPECT_EQ(d_off[b], 4) << "where b = " << b;
    }
  }

  cmprts->need_reorder = true;

  // bnd_post doesn't do the actually final reordering, but
  // let's do it here for a final check
  cmprts->reorder();
  // for (int n = 0; n < cmprts->n_prts; n++) {
  //   float4 xi4 = cmprts->d_xi4[n];
  //   printf("n %d: %g:%g kind %d\n", n, xi4.y, xi4.z, cuda_float_as_int(xi4.w));
  // }
  EXPECT_TRUE(cmprts->check_ordered());

#if 0
  cmprts->dump();
#endif
}

// ======================================================================
// main

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
