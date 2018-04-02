
#include <gtest/gtest.h>

#include "grid.hxx"
#include "fields3d.hxx"
#include "psc_fields_single.h"

struct GridDomain { Grid_t grid; mrc_domain* domain; };

static GridDomain make_grid()
{
  auto domain = Grid_t::Domain{{1, 4, 4},
			       {10.,  40., 40.}, {},
			       {1, 2, 1}};

  mrc_domain* mrc_domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(mrc_domain, "multi");
  mrc_domain_set_param_int3(mrc_domain, "m", domain.gdims);
  mrc_domain_set_param_int3(mrc_domain, "np", domain.np);
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

TEST(Bnd, MakeGrid)
{
  auto grid_domain = make_grid();
  Grid_t& grid = grid_domain.grid;
    
  EXPECT_EQ(grid.domain.gdims, Int3({1, 4, 4}));
  EXPECT_EQ(grid.ldims, Int3({1, 2, 4 }));
  EXPECT_EQ(grid.domain.dx, Grid_t::Real3({ 10., 10., 10. }));
  EXPECT_EQ(grid.n_patches(), 2);
  EXPECT_EQ(grid.patches[0].off, Int3({ 0, 0, 0 }));
  EXPECT_EQ(grid.patches[0].xb, Grid_t::Real3({  0.,  0.,  0. }));
  EXPECT_EQ(grid.patches[0].xe, Grid_t::Real3({ 10., 20., 40. }));
  EXPECT_EQ(grid.patches[1].off, Int3({ 0, 2, 0 }));
  EXPECT_EQ(grid.patches[1].xb, Grid_t::Real3({  0., 20.,  0. }));
  EXPECT_EQ(grid.patches[1].xe, Grid_t::Real3({ 10., 40., 40. }));
}

using Mfields_t = Mfields<fields_single_t>;
const int B = 1;

static void Mfields_dump(Mfields_t& mflds)
{
  for (int p = 0; p < mflds.n_patches(); p++) {
    mflds.grid().Foreach_3d(B, B, [&](int i, int j, int k) {
	for (int m = 0; m < mflds.n_comps(); m++) {
	  printf("p %d ijk [%d:%d:%d] m %d value %02g\n", p, i,j,k, m, mflds[p](m, i,j,k));
	}
      });
  }
}

TEST(Bnd, MakeField)
{
  auto grid_domain = make_grid();
  Grid_t& grid = grid_domain.grid;

  Mfields_t mflds(grid, 1, Int3{ 0, B, B });

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());

  for (int p = 0; p < mflds.n_patches(); p++) {
    int j0 = grid.patches[p].off[1];
    int k0 = grid.patches[p].off[2];
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	int jj = j + j0, kk = k + k0;
	mflds[p](0, i,j,k) = 10*jj + kk;
      });
  }

  Mfields_dump(mflds);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
