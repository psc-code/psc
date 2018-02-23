
#include <gtest/gtest.h>

#include "grid.hxx"

TEST(Grid, Domain)
{
  Int3 gdims = { 8, 4, 2 };
  Int3 ldims = { 4, 2, 2 };
  Grid_t::Real3 length = {  80.,  40., 20. };
  Grid_t::Real3 corner = { -40., -20., 0. };
  std::vector<Int3> offs = { { 0, 0, 0 }, { 4, 0, 0 } };
  Grid_t grid(gdims, ldims, length, corner, offs);

  EXPECT_EQ(grid.gdims, gdims);
  EXPECT_EQ(grid.ldims, ldims);
  EXPECT_EQ(grid.dx, Grid_t::Real3({ 10., 10., 10. }));
  EXPECT_EQ(grid.n_patches(), 2);
  EXPECT_EQ(grid.patches[0].off, offs[0]);
  EXPECT_EQ(grid.patches[0].xb, Grid_t::Real3({ -40., -20.,  0. }));
  EXPECT_EQ(grid.patches[0].xe, Grid_t::Real3({   0.,   0., 20. }));
  EXPECT_EQ(grid.patches[1].off, offs[1]);
  EXPECT_EQ(grid.patches[1].xb, Grid_t::Real3({   0., -20.,  0. }));
  EXPECT_EQ(grid.patches[1].xe, Grid_t::Real3({  40.,   0., 20. }));
}

