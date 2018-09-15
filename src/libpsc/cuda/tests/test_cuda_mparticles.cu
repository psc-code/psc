
#include "cuda_mparticles.h"
#include "cuda_mparticles_sort.cuh"
#include "cuda_collision.cuh"
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

using CudaMparticles = cuda_mparticles<BS144>;

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
    using Real3 = cuda_mparticles_prt::Real3;
    auto ldims = grid_.ldims;
    auto dx = grid_.domain.dx;
    
    int k = n % ldims[2];
    n /= ldims[2];
    int j = n % ldims[1];
    n /= ldims[1];
    int i = n;

    auto prt = cuda_mparticles_prt{
      Real3{Vec3<double>{dx[0] * (i + .5f), dx[1] * (j + .5f), dx[2] * (k + .5f)}},
      Real3{Vec3<int>{i, j, k}},
      1., 0};
    return prt;
  }

private:
  const Grid_t& grid_;
};

void cuda_mparticles_add_particles_test_1(CudaMparticles* cmprts, uint *n_prts_by_patch)
{
  const Grid_t& grid = cmprts->grid_;
  Int3 ldims = grid.ldims;

  for (int p = 0; p < cmprts->n_patches; p++) {
    n_prts_by_patch[p] = ldims[0] * ldims[1] * ldims[2];
  }

  SetParticleTest1 set_particles(grid);
  
  cmprts->reserve_all(n_prts_by_patch);
  cmprts->resize_all(n_prts_by_patch);
  
  for (int p = 0; p < grid.n_patches(); p++) {
    cmprts->set_particles(p, set_particles);
  }
}

// ======================================================================
// CudaMparticlesTest

struct CudaMparticlesTest : TestBase<CudaMparticles>, ::testing::Test
{
  std::unique_ptr<Grid_t> grid_;

  void SetUp()
  {
    auto domain = Grid_t::Domain{{1, 8, 4}, {1., 80., 40.}};
    auto bc = GridBc{};
    auto kinds = Grid_t::Kinds{};
    auto norm = Grid_t::Normalization{};
    double dt = .1;
    grid_.reset(new Grid_t{domain, bc, kinds, norm, dt});
  }
};

// ----------------------------------------------------------------------
TEST_F(CudaMparticlesTest, ConstructorDestructor)
{
  grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
  grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));
  std::unique_ptr<CudaMparticles> cmprts(make_cmprts(*grid_));
  EXPECT_EQ(cmprts->n_patches, 1);
}

// ----------------------------------------------------------------------
TEST_F(CudaMparticlesTest, SetParticles)
{
  grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
  grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));
  std::unique_ptr<CudaMparticles> cmprts(make_cmprts(*grid_));

  uint n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_add_particles_test_1(cmprts.get(), n_prts_by_patch);

  // check that particles are in C order
  int n = 0;
  for (auto prt: cmprts->get_particles(0)) {
    int nn = n++;
    int k = nn % grid_->ldims[2]; nn /= grid_->ldims[2];
    int j = nn % grid_->ldims[1]; nn /= grid_->ldims[1];
    int i = nn;
    EXPECT_FLOAT_EQ(prt.x[0], (i + .5) * grid_->domain.dx[0]);
    EXPECT_FLOAT_EQ(prt.x[1], (j + .5) * grid_->domain.dx[1]);
    EXPECT_FLOAT_EQ(prt.x[2], (k + .5) * grid_->domain.dx[2]);
  }
}

// ---------------------------------------------------------------------
// SetupInternalsDetail
//
// Tests the pieces that go into setup_internals()

TEST_F(CudaMparticlesTest, SetupInternalsDetail)
{
  grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
  grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));

  std::vector<cuda_mparticles_prt> prts = {
    {{ .5, 75., 15. }, {}, 0., 0},
    {{ .5, 35., 15. }, {}, 0., 0},
    {{ .5,  5.,  5. }, {}, 0., 0},
  };
  uint n_prts_by_patch[1];
  n_prts_by_patch[0] = prts.size();

  // can't use make_cmprts() from vector here, since that'll sort etc
  std::unique_ptr<CudaMparticles> cmprts(make_cmprts(*grid_));
  cmprts->reserve_all(n_prts_by_patch);
  cmprts->resize_all(n_prts_by_patch);
  cmprts->set_particles(0, [&](int n) {
      return prts[n];
    });

  auto& d_id = cmprts->by_block_.d_id;
  auto& d_bidx = cmprts->by_block_.d_idx;
  EXPECT_EQ(d_bidx[0], 0);
  EXPECT_EQ(d_bidx[1], 0);
  EXPECT_EQ(d_bidx[2], 0);
  EXPECT_EQ(d_id[0], 0);
  EXPECT_EQ(d_id[1], 0);
  EXPECT_EQ(d_id[2], 0);
  
  EXPECT_TRUE(cmprts->check_in_patch_unordered_slow());
  cmprts->by_block_.find_indices_ids(*cmprts);
  
  EXPECT_EQ(d_bidx[0], 1);
  EXPECT_EQ(d_bidx[1], 0);
  EXPECT_EQ(d_bidx[2], 0);
  EXPECT_EQ(d_id[0], 0);
  EXPECT_EQ(d_id[1], 1);
  EXPECT_EQ(d_id[2], 2);

  EXPECT_TRUE(cmprts->check_bidx_id_unordered_slow());
  cmprts->by_block_.stable_sort();
  
  EXPECT_EQ(d_bidx[0], 0);
  EXPECT_EQ(d_bidx[1], 0);
  EXPECT_EQ(d_bidx[2], 1);
  EXPECT_EQ(d_id[0], 1);
  EXPECT_EQ(d_id[1], 2);
  EXPECT_EQ(d_id[2], 0);

  cmprts->by_block_.reorder_and_offsets(*cmprts);

  float4 xi4_0 = cmprts->d_xi4[0], xi4_1 = cmprts->d_xi4[1], xi4_2 = cmprts->d_xi4[2];
  EXPECT_FLOAT_EQ(xi4_0.y, 35.);
  EXPECT_FLOAT_EQ(xi4_0.z, 15.);
  EXPECT_FLOAT_EQ(xi4_1.y, 5.);
  EXPECT_FLOAT_EQ(xi4_1.z, 5.);
  EXPECT_FLOAT_EQ(xi4_2.y, 75.);
  EXPECT_FLOAT_EQ(xi4_2.z, 15.);

  auto& d_off = cmprts->by_block_.d_off;
  EXPECT_EQ(d_off[0], 0);
  EXPECT_EQ(d_off[1], 2);
  EXPECT_EQ(d_off[2], 3);

  EXPECT_TRUE(cmprts->check_ordered());
}

// ---------------------------------------------------------------------
// SortByCellDetail
//
// Tests the pieces that go into setup_internals()

TEST_F(CudaMparticlesTest, SortByCellDetail)
{
  grid_->kinds.push_back(Grid_t::Kind(-1.,  1., "electron"));
  grid_->kinds.push_back(Grid_t::Kind( 1., 25., "ion"));

  std::vector<cuda_mparticles_prt> prts = {
    {{ .5, 75., 15. }, {}, 0., 0},
    {{ .5, 35., 15. }, {}, 0., 0},
    {{ .5,  5.,  5. }, {}, 0., 0},
  };
  uint n_prts_by_patch[1];
  n_prts_by_patch[0] = prts.size();

  // can't use make_cmprts() from vector here, since that'll sort etc
  std::unique_ptr<CudaMparticles> cmprts(make_cmprts(*grid_));
  cmprts->reserve_all(n_prts_by_patch);
  cmprts->resize_all(n_prts_by_patch);
  cmprts->set_particles(0, [&](int n) {
      return prts[n];
    });
  EXPECT_TRUE(cmprts->check_in_patch_unordered_slow());

  auto sort_by_cell = cuda_mparticles_sort{cmprts->n_cells()};
  auto& d_idx = sort_by_cell.d_idx;
  auto& d_id = sort_by_cell.d_id;
  
  sort_by_cell.find_indices_ids(*cmprts);
  EXPECT_EQ(d_idx[0], 15);
  EXPECT_EQ(d_idx[1], 11);
  EXPECT_EQ(d_idx[2], 0);
  EXPECT_EQ(d_id[0], 0);
  EXPECT_EQ(d_id[1], 1);
  EXPECT_EQ(d_id[2], 2);

  sort_by_cell.stable_sort_cidx();
  EXPECT_EQ(d_idx[0], 0);
  EXPECT_EQ(d_idx[1], 11);
  EXPECT_EQ(d_idx[2], 15);
  EXPECT_EQ(d_id[0], 2);
  EXPECT_EQ(d_id[1], 1);
  EXPECT_EQ(d_id[2], 0);

  sort_by_cell.find_offsets();
  auto& d_off = sort_by_cell.d_off;
  EXPECT_EQ(d_off[0], 0);
  for (int c = 1; c <= 11; c++) {
    EXPECT_EQ(d_off[c], 1) << "c = " << c;
  }
  for (int c = 12; c <= 15; c++) {
    EXPECT_EQ(d_off[c], 2) << "c = " << c;
  }
  EXPECT_EQ(d_off[16], 3);
  
  sort_by_cell.reorder(*cmprts);
  float4 xi4_0 = cmprts->d_xi4[0], xi4_1 = cmprts->d_xi4[1], xi4_2 = cmprts->d_xi4[2];
  EXPECT_FLOAT_EQ(xi4_0.y, 5.);
  EXPECT_FLOAT_EQ(xi4_0.z, 5.);
  EXPECT_FLOAT_EQ(xi4_1.y, 35.);
  EXPECT_FLOAT_EQ(xi4_1.z, 15.);
  EXPECT_FLOAT_EQ(xi4_2.y, 75.);
  EXPECT_FLOAT_EQ(xi4_2.z, 15.);
}

// ----------------------------------------------------------------------
// SetupInternals
//
// tests setup_internals() itself, on a slightly bigger set of particles

TEST_F(CudaMparticlesTest, SetupInternals)
{
  grid_->kinds.push_back(Grid_t::Kind( 1.,  1., "test species"));
  std::unique_ptr<CudaMparticles> cmprts(make_cmprts(*grid_));

  uint n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_add_particles_test_1(cmprts.get(), n_prts_by_patch);

  cmprts->check_in_patch_unordered_slow();

  cmprts->setup_internals();
  
  // check that particles are now in Fortran order
  int cur_bidx = 0;
  for (auto prt: cmprts->get_particles(0)) {
    float4 xi = { prt.x[0], prt.x[1], prt.x[2] };
    int bidx = cmprts->blockIndex(xi, 0);
    EXPECT_GE(bidx, cur_bidx);
    cur_bidx = bidx;
  }

  cmprts->check_ordered();
}

// ----------------------------------------------------------------------
// CudaCollision

TEST_F(CudaMparticlesTest, CudaCollision)
{
  grid_->kinds.push_back(Grid_t::Kind( 1.,  1., "test species"));

  std::vector<cuda_mparticles_prt> prts = {
    {{.5, 75., 15.}, {1.0, 0., 0.}, 1., 0},
    {{.5, 75., 15.}, {1.1, 0., 0.}, 1., 0},
    {{.5, 75., 15.}, {1.2, 0., 0.}, 1., 0},
    {{.5, 35.,  5.}, {0., 1.0, 0.}, 1., 0},
    {{.5, 35.,  5.}, {0., 1.1, 0.}, 1., 0},
    {{.5, 35.,  5.}, {0., 1.2, 0.}, 1., 0},
    {{.5, 35.,  5.}, {0., 1.3, 0.}, 1., 0},
    {{.5, 35.,  5.}, {0., 1.4, 0.}, 1., 0},
    {{.5,  5.,  5.}, {0., 0., 1.0}, 1., 0},
    {{.5,  5.,  5.}, {0., 0., 1.1}, 1., 0},
  };
  uint n_prts_by_patch[1];
  n_prts_by_patch[0] = prts.size();

  // FIXME: can use make_cmprts() from vector here
  std::unique_ptr<CudaMparticles> cmprts(make_cmprts(*grid_));
  cmprts->reserve_all(n_prts_by_patch);
  cmprts->resize_all(n_prts_by_patch);
  cmprts->set_particles(0, [&](int n) {
      return prts[n];
    });

  cmprts->setup_internals();
  cmprts->check_ordered();

  int interval = 1;
  double nu = .1;
  int nicell = 10;
  double dt = .1;
  CudaCollision<CudaMparticles, RngStateCuda> coll{interval, nu, nicell, dt};

  coll(*cmprts);
  for (auto prt: cmprts->get_particles(0)) {
    printf("xi %g %g pxi %g %g %g\n", prt.x[1], prt.x[2],
	   prt.p[0], prt.p[1], prt.p[2]);
  }
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
