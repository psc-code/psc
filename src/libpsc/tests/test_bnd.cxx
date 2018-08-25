
#include <gtest/gtest.h>

#include "grid.hxx"
#include "fields3d.hxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#ifdef USE_CUDA
#include "../libpsc/cuda/bnd_cuda_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_2_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_3_impl.hxx"
#endif

#include "psc_fields_single.h"
#include "psc_fields_c.h"

struct GridDomain { Grid_t grid; mrc_domain* domain; };

static GridDomain make_grid(Int3 gdims, Vec3<double> length)
{
  auto domain = Grid_t::Domain{gdims, length, {},
			       {1, 2, 1}};

  mrc_domain* mrc_domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(mrc_domain, "multi");
  mrc_domain_set_param_int3(mrc_domain, "m", domain.gdims);
  mrc_domain_set_param_int3(mrc_domain, "np", domain.np);
  mrc_domain_set_param_int(mrc_domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(mrc_domain, "bcy", BC_PERIODIC);
  mrc_domain_set_param_int(mrc_domain, "bcz", BC_PERIODIC);

  struct mrc_crds *crds = mrc_domain_get_crds(mrc_domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", 2);
  mrc_crds_set_param_double3(crds, "l", domain.corner);
  mrc_crds_set_param_double3(crds, "h", domain.corner + domain.length);
  
  mrc_domain_setup(mrc_domain);
  //mrc_domain_view(mrc_domain);

  int n_patches;
  auto *mrc_patches = mrc_domain_get_patches(mrc_domain, &n_patches);

  std::vector<Int3> offs;
  for (int p = 0; p < n_patches; p++) {
    offs.emplace_back(mrc_patches[p].off);
  }
  
  return { Grid_t(domain, offs), mrc_domain };
}

template<typename DIM>
static GridDomain make_grid()
{
  return make_grid({1, 8, 4}, {10., 80., 40.});
}

template<>
GridDomain make_grid<dim_xyz>()
{
  return make_grid({2, 8, 4}, {20., 80., 40.});
}

template<typename Mfields>
static void Mfields_dump(Mfields& mflds, int B)
{
  using real_t = typename Mfields::real_t;
  
  for (int p = 0; p < mflds.n_patches(); p++) {
    mflds.grid().Foreach_3d(B, B, [&](int i, int j, int k) {
	for (int m = 0; m < mflds.n_comps(); m++) {
	  printf("p %d ijk [%d:%d:%d] m %d value %02g\n", p, i,j,k, m, (real_t) mflds[p](m, i,j,k));
	}
      });
  }
}

TEST(Bnd, MakeGrid)
{
  auto grid_domain = make_grid<dim_yz>();
  Grid_t& grid = grid_domain.grid;
    
  EXPECT_EQ(grid.domain.gdims, Int3({1, 8, 4}));
  EXPECT_EQ(grid.ldims, Int3({1, 4, 4 }));
  EXPECT_EQ(grid.domain.dx, Grid_t::Real3({ 10., 10., 10. }));

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_size == 1) {
    EXPECT_EQ(grid.n_patches(), 2);
    EXPECT_EQ(grid.patches[0].off, Int3({ 0, 0, 0 }));
    EXPECT_EQ(grid.patches[0].xb, Grid_t::Real3({  0.,  0.,  0. }));
    EXPECT_EQ(grid.patches[0].xe, Grid_t::Real3({ 10., 40., 40. }));
    EXPECT_EQ(grid.patches[1].off, Int3({ 0, 4, 0 }));
    EXPECT_EQ(grid.patches[1].xb, Grid_t::Real3({  0., 40.,  0. }));
    EXPECT_EQ(grid.patches[1].xe, Grid_t::Real3({ 10., 80., 40. }));
  }
}

template<typename BND, typename DIM>
struct TestConfigBnd
{
  using Bnd = BND;
  using dim = DIM;
};

template<typename T>
struct BndTest : public ::testing::Test
{
  using Bnd = typename T::Bnd;
  using dim = typename T::dim;
};

using BndTestTypes = ::testing::Types<TestConfigBnd<Bnd_<MfieldsSingle>, dim_yz>,
				      TestConfigBnd<Bnd_<MfieldsC>, dim_yz>,
#ifdef USE_CUDA
				      TestConfigBnd<BndCuda, dim_yz>,
				      TestConfigBnd<BndCuda2<MfieldsCuda>, dim_xyz>,
				      TestConfigBnd<BndCuda3<MfieldsCuda>, dim_xyz>,
#endif
				      TestConfigBnd<Bnd_<MfieldsSingle>, dim_xyz>>;

TYPED_TEST_CASE(BndTest, BndTestTypes);

const int B = 2;

TYPED_TEST(BndTest, FillGhosts)
{
  using Base = BndTest<TypeParam>;
  using Bnd = typename Base::Bnd;
  using dim = typename Base::dim;
  using Mfields = typename Bnd::Mfields;

  auto grid_domain = make_grid<dim>();
  auto& grid = grid_domain.grid;
  auto ibn = Int3{B, B, B};
  if (dim::InvarX::value) ibn[0] = 0;
  auto mflds = Mfields{grid, 1, ibn};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  for (int p = 0; p < mflds.n_patches(); p++) {
    int i0 = grid.patches[p].off[0];
    int j0 = grid.patches[p].off[1];
    int k0 = grid.patches[p].off[2];
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	int ii = i + i0, jj = j + j0, kk = k + k0;
	mflds[p](0, i,j,k) = 100*ii + 10*jj + kk;
      });
  }

  //Mfields_dump(mflds, B);
  for (int p = 0; p < mflds.n_patches(); p++) {
    int i0 = grid.patches[p].off[0];
    int j0 = grid.patches[p].off[1];
    int k0 = grid.patches[p].off[2];
    grid.Foreach_3d(B, B, [&](int i, int j, int k) {
	int ii = i + i0, jj = j + j0, kk = k + k0;
	if (i >= 0 && i < grid.ldims[0] &&
	    j >= 0 && j < grid.ldims[1] &&
	    k >= 0 && k < grid.ldims[2]) {
	  EXPECT_EQ(mflds[p](0, i,j,k), 100*ii + 10*jj + kk);
	} else {
	  EXPECT_EQ(mflds[p](0, i,j,k), 0);
	}
      });
  }

  Bnd bnd{grid, ibn};
  bnd.fill_ghosts(mflds, 0, 1);

  //Mfields_dump(mflds, B);
  for (int p = 0; p < mflds.n_patches(); p++) {
    int i0 = grid.patches[p].off[0];
    int j0 = grid.patches[p].off[1];
    int k0 = grid.patches[p].off[2];
    grid.Foreach_3d(B, B, [&](int i, int j, int k) {
	int ii = i + i0, jj = j + j0, kk = k + k0;
	ii = (ii + grid.domain.gdims[0]) % grid.domain.gdims[0];
	jj = (jj + grid.domain.gdims[1]) % grid.domain.gdims[1];
	kk = (kk + grid.domain.gdims[2]) % grid.domain.gdims[2];
	EXPECT_EQ(mflds[p](0, i,j,k), 100*ii + 10*jj + kk);
      });
  }

  // let's do it again to test CudaBnd caching
  bnd.fill_ghosts(mflds, 0, 1);
}

TYPED_TEST(BndTest, AddGhosts)
{
  using Base = BndTest<TypeParam>;
  using Bnd = typename Base::Bnd;
  using dim = typename Base::dim;
  using Mfields = typename Bnd::Mfields;

  auto grid_domain = make_grid<dim>();
  auto& grid = grid_domain.grid;
  auto ibn = Int3{B, B, B};
  if (dim::InvarX::value) ibn[0] = 0;
  auto mflds = Mfields{grid, 1, ibn};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  for (int p = 0; p < mflds.n_patches(); p++) {
    int i0 = grid.patches[p].off[0];
    int j0 = grid.patches[p].off[1];
    int k0 = grid.patches[p].off[2];
    grid.Foreach_3d(B, B, [&](int i, int j, int k) {
	mflds[p](0, i,j,k) = 1;
      });
  }

  //Mfields_dump(mflds, B);

  Bnd bnd{grid, ibn};
  bnd.add_ghosts(mflds, 0, 1);

  //Mfields_dump(mflds, 0*B);
  for (int p = 0; p < mflds.n_patches(); p++) {
    int j0 = grid.patches[p].off[1];
    int k0 = grid.patches[p].off[2];
    auto& ldims = grid.ldims;
    grid.Foreach_3d(B, B, [&](int i, int j, int k) {
	int n_nei = 0;
	if (i >= 0 && i < ldims[0] &&
	    j >= 0 && j < ldims[1] &&
	    k >= 0 && k < ldims[2]) {
	  if (!dim::InvarX::value) {
	    if (i < B || i >= ldims[0] - B) n_nei++;
	  }
	  if (j < B || j >= ldims[1] - B) n_nei++;
	  if (k < B || k >= ldims[2] - B) n_nei++;
	  if (n_nei == 2) n_nei = 3; // edge -> diagonal contribution, too
	  else if (n_nei == 3) n_nei = 11; // corner -> why? FIXME why not 7?
	}
  	EXPECT_EQ(mflds[p](0, i,j,k), 1 + n_nei) << "ijk " << i << " " << j << " " << k;
      });
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
