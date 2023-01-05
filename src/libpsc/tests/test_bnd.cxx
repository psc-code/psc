
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

static Grid_t make_grid(Int3 gdims, Vec3<double> length)
{
  auto domain = Grid_t::Domain{gdims, length, {}, {1, 2, 1}};
  auto bc = psc::grid::BC{};
  auto kinds = Grid_t::Kinds{};
  auto norm = Grid_t::Normalization{};
  double dt = .1;
  int n_patches = -1;

  return Grid_t{domain, bc, kinds, norm, dt, n_patches, {B, B, B}};
}

template <typename DIM>
static Grid_t make_grid()
{
  return make_grid({1, 8, 4}, {10., 80., 40.});
}

template <>
Grid_t make_grid<dim_xyz>()
{
  return make_grid({2, 8, 4}, {20., 80., 40.});
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
  ::testing::Types<TestConfigBnd<MfieldsSingle, Bnd_<MfieldsSingle>, dim_yz>,
                   TestConfigBnd<MfieldsC, Bnd_<MfieldsC>, dim_yz>,
#ifdef USE_CUDA
                   TestConfigBnd<MfieldsCuda, BndCuda3<MfieldsCuda>, dim_xyz>,
                   TestConfigBnd<MfieldsCuda, Bnd_<MfieldsCuda>, dim_xyz>,
#endif
                   TestConfigBnd<MfieldsSingle, Bnd_<MfieldsSingle>, dim_xyz>>;

TYPED_TEST_SUITE(BndTest, BndTestTypes);

TYPED_TEST(BndTest, FillGhosts)
{
  using Mfields = typename TypeParam::Mfields;
  using Bnd = typename TypeParam::Bnd;
  using dim = typename TypeParam::dim;

  auto grid = make_grid<dim>();
  auto ibn = Int3{B, B, B};
  if (dim::InvarX::value)
    ibn[0] = 0;
  auto mflds = Mfields{grid, 1, ibn};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  {
    auto&& h_mflds = gt::host_mirror(mflds.gt());
    h_mflds.view() = 0.;
    for (int p = 0; p < mflds.n_patches(); p++) {
      int i0 = grid.patches[p].off[0];
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      auto flds =
        make_Fields3d<dim_xyz>(h_mflds.view(_all, _all, _all, _all, p), -ibn);
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
        int ii = i + i0, jj = j + j0, kk = k + k0;
        flds(0, i, j, k) = 100 * ii + 10 * jj + kk;
      });
    }
    gt::copy(h_mflds, mflds.storage());
  }

  // Mfields_dump(mflds, B);
  {
    auto&& h_mflds = gt::host_mirror(mflds.gt());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds =
        make_Fields3d<dim_xyz>(h_mflds.view(_all, _all, _all, _all, p), -ibn);
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

  Bnd bnd{grid, ibn};
  bnd.fill_ghosts(mflds, 0, 1);

  // Mfields_dump(mflds, B);
  {
    auto&& h_mflds = gt::host_mirror(mflds.gt());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds =
        make_Fields3d<dim_xyz>(h_mflds.view(_all, _all, _all, _all, p), -ibn);
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
  auto ibn = Int3{B, B, B};
  if (dim::InvarX::value)
    ibn[0] = 0;
  auto mflds = Mfields{grid, 1, ibn};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  {
    auto&& h_mflds = gt::host_mirror(mflds.gt());
    h_mflds.view() = 0.;
    for (int p = 0; p < mflds.n_patches(); p++) {
      int i0 = grid.patches[p].off[0];
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      auto flds =
        make_Fields3d<dim_xyz>(h_mflds.view(_all, _all, _all, _all, p), -ibn);
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
        int ii = i + i0, jj = j + j0, kk = k + k0;
        flds(0, i, j, k) = 100 * ii + 10 * jj + kk;
      });
    }
    gt::copy(h_mflds, mflds.storage());
  }

  Bnd bnd{grid, ibn};
  bnd.fill_ghosts(mflds.grid(), mflds.storage(), mflds.ib(), 0, 1);

  {
    auto&& h_mflds = gt::host_mirror(mflds.gt());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds =
        make_Fields3d<dim_xyz>(h_mflds.view(_all, _all, _all, _all, p), -ibn);
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
  auto ibn = Int3{B, B, B};
  if (dim::InvarX::value)
    ibn[0] = 0;
  auto mflds = Mfields{grid, 1, ibn};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  {
    auto&& h_mflds = gt::host_mirror(mflds.gt());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds =
        make_Fields3d<dim_xyz>(h_mflds.view(_all, _all, _all, _all, p), -ibn);
      int i0 = grid.patches[p].off[0];
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      grid.Foreach_3d(B, B, [&](int i, int j, int k) { flds(0, i, j, k) = 1; });
    }
    gt::copy(h_mflds, mflds.storage());
  }

  // Mfields_dump(mflds, B);

  Bnd bnd{grid, ibn};
  bnd.add_ghosts(mflds, 0, 1);

  // Mfields_dump(mflds, 0*B);
  {
    auto&& h_mflds = gt::host_mirror(mflds.gt());
    gt::copy(mflds.storage(), h_mflds);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto flds =
        make_Fields3d<dim_xyz>(h_mflds.view(_all, _all, _all, _all, p), -ibn);
      int j0 = grid.patches[p].off[1];
      int k0 = grid.patches[p].off[2];
      auto& ldims = grid.ldims;
      grid.Foreach_3d(B, B, [&](int i, int j, int k) {
        int n_nei = 0;
        if (i >= 0 && i < ldims[0] && j >= 0 && j < ldims[1] && k >= 0 &&
            k < ldims[2]) {
          if (!dim::InvarX::value) {
            if (i < B || i >= ldims[0] - B)
              n_nei++;
          }
          if (j < B || j >= ldims[1] - B)
            n_nei++;
          if (k < B || k >= ldims[2] - B)
            n_nei++;
          if (n_nei == 2)
            n_nei = 3; // edge -> diagonal contribution, too
          else if (n_nei == 3)
            n_nei = 11; // corner -> why? FIXME why not 7?
        }
        EXPECT_EQ(flds(0, i, j, k), 1 + n_nei)
          << "ijk " << i << " " << j << " " << k;
      });
    }
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
