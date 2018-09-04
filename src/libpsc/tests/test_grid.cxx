
#include <gtest/gtest.h>

#include "test_common.hxx"

TEST(Grid, Domain)
{
  Grid_t grid = MakeTestGrid{}();
    
  EXPECT_EQ(grid.domain.gdims, Int3({ 8, 4, 2 }));
  EXPECT_EQ(grid.ldims, Int3({ 4, 2, 2 }));
  EXPECT_EQ(grid.domain.dx, Grid_t::Real3({ 10., 10., 10. }));
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
  Grid_t grid = MakeTestGrid{}();

  grid.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  EXPECT_EQ(grid.kinds.size(), 1);
}

