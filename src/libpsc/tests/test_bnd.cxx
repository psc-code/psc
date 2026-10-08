
#include <gtest/gtest.h>

#include "grid.hxx"
#include "fields3d.hxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#ifdef USE_CUDA
#include "../libpsc/cuda/bnd_cuda_3_impl.hxx"
#endif

#include "psc_fields_single.h"
#include "psc_fields_c.h"

const int B = 2;

template <typename D>
static Grid_t make_grid(Int3 gdims, Vec3<double> length)
{
  auto domain = Grid_t::Domain{gdims, length, {}, {1, 2, 1}};
  auto bc = psc::grid::BC{};
  auto kinds = Grid_t::Kinds{};
  auto norm = Grid_t::Normalization{};
  double dt = .1;
  int n_patches = -1;

  auto ibn = Int3{B, B, B};
  if (D::InvarX::value) {
    ibn[0] = 0;
  }
  return Grid_t{domain, bc, kinds, norm, dt, n_patches, ibn};
}

template <typename D>
static Grid_t make_grid()
{
  return make_grid<D>({1, 8, 4}, {10., 80., 40.});
}

template <>
Grid_t make_grid<dim_xyz>()
{
  return make_grid<dim_xyz>({2, 8, 4}, {20., 80., 40.});
}

template <typename Mfields>
static void Mfields_dump(Mfields& mflds, int B)
{
  using real_t = typename Mfields::real_t;

  for (int p = 0; p < mflds.n_patches(); p++) {
    mflds.grid().Foreach_3d(B, B, [&](int i, int j, int k) {
      for (int m = 0; m < mflds.n_comps(); m++) {
        printf("p %d ijk [%d:%d:%d] m %d value %02g\n", p, i, j, k, m,
               (real_t)mflds[p](m, i, j, k));
      }
    });
  }
}

TEST(Bnd, MakeGrid)
{
  auto grid = make_grid<dim_yz>();

  EXPECT_EQ(grid.domain.gdims, Int3({1, 8, 4}));
  EXPECT_EQ(grid.ldims, Int3({1, 4, 4}));
  EXPECT_EQ(grid.domain.dx, Grid_t::Real3({10., 10., 10.}));

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_size == 1) {
    EXPECT_EQ(grid.n_patches(), 2);
    EXPECT_EQ(grid.patches[0].off, Int3({0, 0, 0}));
    EXPECT_EQ(grid.patches[0].xb, Grid_t::Real3({0., 0., 0.}));
    EXPECT_EQ(grid.patches[0].xe, Grid_t::Real3({10., 40., 40.}));
    EXPECT_EQ(grid.patches[1].off, Int3({0, 4, 0}));
    EXPECT_EQ(grid.patches[1].xb, Grid_t::Real3({0., 40., 0.}));
    EXPECT_EQ(grid.patches[1].xe, Grid_t::Real3({10., 80., 40.}));
  }
}

template <typename MF, typename BND, typename DIM>
struct TestConfigBnd
{
  using Mfields = MF;
  using Bnd = BND;
  using dim = DIM;
};

template <typename T>
struct BndTest : public ::testing::Test
{};

using BndTestTypes =
  ::testing::Types<TestConfigBnd<MfieldsSingle, Bnd_, dim_yz>,
                   TestConfigBnd<MfieldsC, Bnd_, dim_yz>,
#ifdef USE_CUDA
                   TestConfigBnd<MfieldsCuda, BndCuda3, dim_xyz>,
                   TestConfigBnd<MfieldsCuda, Bnd_, dim_xyz>,
#endif
                   TestConfigBnd<MfieldsSingle, Bnd_, dim_xyz>>;

TYPED_TEST_SUITE(BndTest, BndTestTypes);

TYPED_TEST(BndTest, FillGhosts)
{
  using Mfields = typename TypeParam::Mfields;
  using Bnd = typename TypeParam::Bnd;
  using dim = typename TypeParam::dim;

  auto grid = make_grid<dim>();
  auto mflds = Mfields{grid, 1, grid.ibn};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  {
    auto&& h_mflds = gt::host_mirror(mflds.storage());
    h_mflds.view() = 0.;
    for (int p = 0; p < mflds.n_patches(); p++) {
      int i0 = grid.patches[p].off[0];
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      auto flds = make_Fields3d<dim_xyz>(
        h_mflds.view(_all, _all, _all, _all, p), -grid.ibn);
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
        int ii = i + i0, jj = j + j0, kk = k + k0;
        flds(0, i, j, k) = 100 * ii + 10 * jj + kk;
      });
    }
    gt::copy(h_mflds, mflds.storage());
  }

  // Mfields_dump(mflds, B);
  {
    auto&& h_mflds = gt::host_mirror(mflds.storage());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds = make_Fields3d<dim_xyz>(
        h_mflds.view(_all, _all, _all, _all, p), -grid.ibn);
      int i0 = grid.patches[p].off[0];
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      grid.Foreach_3d(B, B, [&](int i, int j, int k) {
        int ii = i + i0, jj = j + j0, kk = k + k0;
        if (i >= 0 && i < grid.ldims[0] && j >= 0 && j < grid.ldims[1] &&
            k >= 0 && k < grid.ldims[2]) {
          EXPECT_EQ(flds(0, i, j, k), 100 * ii + 10 * jj + kk);
        } else {
          EXPECT_EQ(flds(0, i, j, k), 0);
        }
      });
    }
  }

  Bnd bnd;
  bnd.fill_ghosts(mflds, 0, 1);

  // Mfields_dump(mflds, B);
  {
    auto&& h_mflds = gt::host_mirror(mflds.storage());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds = make_Fields3d<dim_xyz>(
        h_mflds.view(_all, _all, _all, _all, p), -grid.ibn);
      int i0 = grid.patches[p].off[0];
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      grid.Foreach_3d(B, B, [&](int i, int j, int k) {
        int ii = i + i0, jj = j + j0, kk = k + k0;
        ii = (ii + grid.domain.gdims[0]) % grid.domain.gdims[0];
        jj = (jj + grid.domain.gdims[1]) % grid.domain.gdims[1];
        kk = (kk + grid.domain.gdims[2]) % grid.domain.gdims[2];
        EXPECT_EQ(flds(0, i, j, k), 100 * ii + 10 * jj + kk);
      });
    }
  }

  // let's do it again to test CudaBnd caching
  bnd.fill_ghosts(mflds, 0, 1);
}

// almost same as "FillGhosts" but uses gt-based interface bnd
TYPED_TEST(BndTest, FillGhostsGt)
{
  using Mfields = typename TypeParam::Mfields;
  using Bnd = typename TypeParam::Bnd;
  using dim = typename TypeParam::dim;

  auto grid = make_grid<dim>();
  auto mflds = Mfields{grid, 1, grid.ibn};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  {
    auto&& h_mflds = gt::host_mirror(mflds.storage());
    h_mflds.view() = 0.;
    for (int p = 0; p < mflds.n_patches(); p++) {
      int i0 = grid.patches[p].off[0];
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      auto flds = make_Fields3d<dim_xyz>(
        h_mflds.view(_all, _all, _all, _all, p), -grid.ibn);
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
        int ii = i + i0, jj = j + j0, kk = k + k0;
        flds(0, i, j, k) = 100 * ii + 10 * jj + kk;
      });
    }
    gt::copy(h_mflds, mflds.storage());
  }

  Bnd bnd{};
  bnd.fill_ghosts(mflds.grid(), mflds.storage(), mflds.ib(), 0, 1);

  {
    auto&& h_mflds = gt::host_mirror(mflds.storage());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds = make_Fields3d<dim_xyz>(
        h_mflds.view(_all, _all, _all, _all, p), -grid.ibn);
      int i0 = grid.patches[p].off[0];
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      grid.Foreach_3d(B, B, [&](int i, int j, int k) {
        int ii = i + i0, jj = j + j0, kk = k + k0;
        ii = (ii + grid.domain.gdims[0]) % grid.domain.gdims[0];
        jj = (jj + grid.domain.gdims[1]) % grid.domain.gdims[1];
        kk = (kk + grid.domain.gdims[2]) % grid.domain.gdims[2];
        EXPECT_EQ(flds(0, i, j, k), 100 * ii + 10 * jj + kk);
      });
    }
  }
}

TYPED_TEST(BndTest, AddGhosts)
{
  using Mfields = typename TypeParam::Mfields;
  using Bnd = typename TypeParam::Bnd;
  using dim = typename TypeParam::dim;

  auto grid = make_grid<dim>();
  auto mflds = Mfields{grid, 1, grid.ibn};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  {
    auto&& h_mflds = gt::host_mirror(mflds.storage());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds = make_Fields3d<dim_xyz>(
        h_mflds.view(_all, _all, _all, _all, p), -grid.ibn);

      grid.Foreach_3d(B, B, [&](int i, int j, int k) { flds(0, i, j, k) = 1; });
    }
    gt::copy(h_mflds, mflds.storage());
  }

  // Mfields_dump(mflds, B);

  Bnd bnd{};
  bnd.add_ghosts(mflds, 0, 1);

  // Mfields_dump(mflds, 0*B);
  {
    auto&& h_mflds = gt::host_mirror(mflds.storage());
    gt::copy(mflds.storage(), h_mflds);

    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds = make_Fields3d<dim_xyz>(
        h_mflds.view(_all, _all, _all, _all, p), -grid.ibn);
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      auto& ldims = grid.ldims;

      grid.Foreach_3d(B, B, [&](int i, int j, int k) {
        int n_neighbors_x = 0;
        int n_neighbors_y = 0;
        int n_neighbors_z = 0;

        if (i >= 0 && i < ldims[0] && j >= 0 && j < ldims[1] && k >= 0 &&
            k < ldims[2]) {
          if (!dim::InvarX::value) {
            n_neighbors_x += i < B;
            n_neighbors_x += i >= ldims[0] - B;
          }

          n_neighbors_y += j < B;
          n_neighbors_y += j >= ldims[1] - B;

          n_neighbors_z += k < B;
          n_neighbors_z += k >= ldims[2] - B;
        }

        int n_neighbors =
          (n_neighbors_x + 1) * (n_neighbors_y + 1) * (n_neighbors_z + 1) - 1;

        EXPECT_EQ(flds(0, i, j, k), 1 + n_neighbors)
          << "ijk " << i << " " << j << " " << k;
      });
    }
  }
}

// ======================================================================
// Non-periodic y, 2 periodic patches in z. Since there are no neighbors in y,
// the y-ghost rows are exchanged along z like interior rows; in particular, the
// (y-ghost, z-ghost) corners must be communicated.

static Grid_t make_grid_open_y()
{
  auto domain = Grid_t::Domain{{1, 4, 16}, {10., 40., 160.}, {}, {1, 1, 2}};
  auto bc = psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_PERIODIC},
                          {BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_PERIODIC},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_PERIODIC},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_PERIODIC}};
  auto kinds = Grid_t::Kinds{};
  auto norm = Grid_t::Normalization{};
  double dt = .1;
  int n_patches = -1;
  auto ibn = Int3{0, B, B};
  return Grid_t{domain, bc, kinds, norm, dt, n_patches, ibn};
}

TEST(Bnd, FillGhostsOpenY)
{
  auto grid = make_grid_open_y();
  auto mflds = MfieldsC{grid, 1, grid.ibn};
  auto& ldims = grid.ldims;
  auto& gdims = grid.domain.gdims;

  // y ghosts are not touched by the exchange along y, so give them values too
  auto value = [&](int jj, int kk) { return 100 * (jj + B) + kk; };

  for (int p = 0; p < mflds.n_patches(); p++) {
    auto flds = make_Fields3d<dim_xyz>(mflds[p]);
    Int3 off = grid.patches[p].off;
    grid.Foreach_3d(B, B, [&](int i, int j, int k) {
      bool interior_z = k >= 0 && k < ldims[2];
      flds(0, i, j, k) = interior_z ? value(j + off[1], k + off[2]) : 0;
    });
  }

  Bnd_ bnd;
  bnd.fill_ghosts(mflds, 0, 1);

  for (int p = 0; p < mflds.n_patches(); p++) {
    auto flds = make_Fields3d<dim_xyz>(mflds[p]);
    Int3 off = grid.patches[p].off;
    grid.Foreach_3d(B, B, [&](int i, int j, int k) {
      int kk = (k + off[2] + gdims[2]) % gdims[2];
      EXPECT_EQ(flds(0, i, j, k), value(j + off[1], kk))
        << "p " << p << " jk " << j << " " << k;
    });
  }
}

TEST(Bnd, AddGhostsOpenY)
{
  auto grid = make_grid_open_y();
  auto mflds = MfieldsC{grid, 1, grid.ibn};
  auto& ldims = grid.ldims;

  for (int p = 0; p < mflds.n_patches(); p++) {
    auto flds = make_Fields3d<dim_xyz>(mflds[p]);
    grid.Foreach_3d(B, B, [&](int i, int j, int k) { flds(0, i, j, k) = 1; });
  }

  Bnd_ bnd;
  bnd.add_ghosts(mflds, 0, 1);

  for (int p = 0; p < mflds.n_patches(); p++) {
    auto flds = make_Fields3d<dim_xyz>(mflds[p]);
    grid.Foreach_3d(B, B, [&](int i, int j, int k) {
      bool interior_z = k >= 0 && k < ldims[2];
      bool near_edge_z = k < B || k >= ldims[2] - B;
      int expected = 1 + (interior_z && near_edge_z);
      EXPECT_EQ(flds(0, i, j, k), expected)
        << "p " << p << " jk " << j << " " << k;
    });
  }
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
