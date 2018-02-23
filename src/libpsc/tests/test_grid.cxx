
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

