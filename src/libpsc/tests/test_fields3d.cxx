
#include <gtest/gtest.h>

#include <fields3d.hxx>


using Real = double;
using Layout = kg::LayoutSOA;

static Grid_t make_grid()
{
  auto domain =
    Grid_t::Domain{{8, 4, 2}, {80., 40., 20.}, {-40., -20., 0.}, {2, 2, 1}};
  auto bc = GridBc{};
  auto kinds = Grid_t::Kinds{};
  auto norm = Grid_t::Normalization{};
  double dt = .1;
  return Grid_t{domain, bc, kinds, norm, dt};
}

TEST(fields3d, CtorNoArray)
{
  auto grid = make_grid();
  auto f = fields3d<Real, Layout>{grid, {1, 2, 3}, {6, 5, 4}, 2};
  f.dtor();
}

TEST(fields3d, CtorArray)
{
  auto grid = make_grid();
  auto storage = std::vector<Real>(5 * 3 * 1 * 2);
  auto f = fields3d<Real, Layout>{grid, {1, 2, 3}, {6, 5, 4}, 2, storage.data()};
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
