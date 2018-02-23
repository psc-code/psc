
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

#include "psc.h" // FIXME, just for EX, ...
#include "psc_fields_single.h"

template<typename mfields, typename Set>
void setValues(mfields& mflds, Set set)
{
  for (int p = 0; p < mflds.n_patches(); p++) {
    auto F = mflds[p];
    F(EX, 0,0,0) = set(EX);
    F(EX, 0,1,0) = set(EX);
    F(EX, 0,0,1) = set(EX);
    F(EX, 0,1,1) = set(EX);

    F(EY, 0,0,0) = set(EY);
    F(EY, 0,0,1) = set(EY);
    F(EY, 1,0,0) = set(EY);
    F(EY, 1,0,1) = set(EY);
    
    F(EZ, 0,0,0) = set(EZ);
    F(EZ, 1,0,0) = set(EZ);
    F(EZ, 0,1,0) = set(EZ);
    F(EZ, 1,1,0) = set(EZ);

    F(HX, 0,0,0) = set(HX);
    F(HX, 1,0,0) = set(HX);

    F(HY, 0,0,0) = set(HY);
    F(HY, 0,1,0) = set(HY);

    F(HZ, 0,0,0) = set(HZ);
    F(HZ, 0,0,1) = set(HZ);
  }
}
  
TEST(mfields, Constructor)
{
  using Mfields = psc_mfields_<fields_single_t>;

  Grid_t grid = make_grid();
  Mfields mflds(grid, NR_FIELDS, Int3{ 1, 1, 1 });

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());
}

TEST(mfields, Set)
{
  using Mfields = psc_mfields_<fields_single_t>;

  Grid_t grid = make_grid();
  Mfields mflds(grid, NR_FIELDS, Int3{ 1, 1, 1 });

  setValues(mflds, [](int m) -> Mfields::real_t {
      switch(m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  auto F = mflds[0];
  EXPECT_EQ(F(EY, 0, 0, 0), 2.);
}

