
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
  auto f = fields3d<Real, Layout>{grid, {1, 2, 3}, {2, 3, 4}, 2};
  f.dtor();
}

TEST(fields3d, CtorArray)
{
  auto grid = make_grid();
  auto storage = std::vector<Real>(2 * 3 * 4 * 2);
  auto f = fields3d<Real, Layout>{grid, {1, 2, 3}, {2, 3, 4}, 2, storage.data()};
}

TEST(fields3d, BoundsEtc)
{
  auto grid = make_grid();
  auto storage = std::vector<Real>(2 * 3 * 4 * 2);
  auto f = fields3d<Real, Layout>{grid, {1, 2, 3}, {2, 3, 4}, 2, storage.data()};

  EXPECT_EQ(f.ib(), Int3({1, 2, 3}));
  EXPECT_EQ(f.im(), Int3({2, 3, 4}));
  EXPECT_EQ(f.n_comps(), 2);
  EXPECT_EQ(f.n_cells(), 2 * 3 * 4);
  EXPECT_EQ(f.size(), 2 * 3 * 4 * 2);
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
