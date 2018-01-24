
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
    grid_->gdims = { 1, 8, 4 };
    grid_->ldims = { 1, 4, 2 };
    grid_->dx = Double3{ 1., 80., 40. } / Double3(grid_->gdims);
    grid_->patches.emplace_back(Grid_t::Patch({ 0.,  0.,  0.}, { 1., 40., 20. }));
    grid_->patches.emplace_back(Grid_t::Patch({ 0., 40.,  0.}, { 1., 80., 20. }));
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

  cmprts->dump();
}