
#include <gtest/gtest.h>

#include <fields3d.hxx>

#include <kg/SArray.h>

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

TEST(SArray, Ctor)
{
  auto f = kg::SArray<Real, Layout>{{1, 2, 3}, {2, 3, 4}, 2};
}

TEST(fields3d_view, Ctor)
{
  auto grid = make_grid();
  auto storage = std::vector<Real>(2 * 3 * 4 * 2);
  auto f = fields3d_view<Real, Layout>{grid, {1, 2, 3}, {2, 3, 4}, 2, storage.data()};
}

TEST(SArray, BoundsEtc)
{
  auto f = kg::SArray<Real, Layout>{{1, 2, 3}, {2, 3, 4}, 2};

  EXPECT_EQ(f.ib(), Int3({1, 2, 3}));
  EXPECT_EQ(f.im(), Int3({2, 3, 4}));
  EXPECT_EQ(f.n_comps(), 2);
  EXPECT_EQ(f.n_cells(), 2 * 3 * 4);
  EXPECT_EQ(f.size(), 2 * 3 * 4 * 2);
}

TEST(SArray, index)
{
  auto f = kg::SArray<Real, Layout>{{1, 2, 3}, {2, 3, 4}, 2};

  EXPECT_EQ(f.index(0, 1, 2, 3), 0);
  EXPECT_EQ(f.index(0, 2, 2, 3), 1);
  EXPECT_EQ(f.index(1, 1, 2, 3), 2 * 3 * 4);
}

TEST(SArray, data)
{
  auto f = kg::SArray<Real, Layout>{{1, 2, 3}, {2, 3, 4}, 2};

  EXPECT_EQ(f.data(), &f(0, 1, 2, 3));
  const auto& fc = const_cast<const kg::SArray<Real, Layout>&>(f);
  EXPECT_EQ(fc.data(), &fc(0, 1, 2, 3));
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
