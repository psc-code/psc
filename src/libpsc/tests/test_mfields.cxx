
#include <gtest/gtest.h>

#include "fields3d.hxx"
#include "psc_fields_c.h"
#include "psc_fields_single.h"
#ifdef USE_CUDA
#include "psc_fields_cuda.h"
#endif
#include "setup_fields.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#endif

#include "psc.h" // FIXME, just for EX etc

static Grid_t make_grid()
{
  auto domain =
    Grid_t::Domain{{8, 4, 2}, {80., 40., 20.}, {-40., -20., 0.}, {2, 2, 1}};
  auto bc = psc::grid::BC{};
  auto kinds = Grid_t::Kinds{};
  auto norm = Grid_t::Normalization{};
  double dt = .1;
  return Grid_t{domain, bc, kinds, norm, dt};
}

template <typename T>
class MfieldsTest : public ::testing::Test
{};

#ifdef xUSE_CUDA
using MfieldsTestTypes = ::testing::Types<MfieldsSingle, MfieldsC, MfieldsCuda>;
#else
using MfieldsTestTypes = ::testing::Types<MfieldsSingle, MfieldsC>;
#endif

TYPED_TEST_SUITE(MfieldsTest, MfieldsTestTypes);

TYPED_TEST(MfieldsTest, Constructor)
{
  using Mfields = TypeParam;

  auto grid = make_grid();
  auto mflds = Mfields{grid, NR_FIELDS, Int3{1, 1, 1}};

  EXPECT_EQ(mflds.n_patches(), grid.n_patches());
}

TYPED_TEST(MfieldsTest, Access)
{
  using Mfields = TypeParam;

  auto grid = make_grid();
  auto mflds = Mfields{grid, NR_FIELDS, Int3{1, 1, 1}};
  auto flds = make_Fields3d<dim_xyz>(mflds[0]);

  EXPECT_EQ(flds(0, 1, 1, 1), 0.);

  flds(0, 1, 1, 1) = 99.;

  EXPECT_EQ(flds(0, 1, 1, 1), 99.);
}

TYPED_TEST(MfieldsTest, ZeroComp)
{
  using Mfields = TypeParam;

  auto grid = make_grid();
  auto mflds = Mfields{grid, NR_FIELDS, Int3{1, 1, 1}};
  auto flds = make_Fields3d<dim_xyz>(mflds[0]);

  flds(EX, 4, 2, 2) = 1.;
  flds(EY, -1, -1, -1) = 2.;
  flds(EY, 0, 0, 0) = 3.;
  flds(EY, 4, 2, 2) = 4.;
  flds(EZ, -1, -1, -1) = 5.;

  mflds.storage().view(_all, _all, _all, EY, _all) = 0.;

  EXPECT_EQ(flds(EX, 4, 2, 2), 1.);
  EXPECT_EQ(flds(EY, -1, -1, -1), 0.);
  EXPECT_EQ(flds(EY, 0, 0, 0), 0.);
  EXPECT_EQ(flds(EY, 4, 2, 2), 0.);
  EXPECT_EQ(flds(EZ, -1, -1, -1), 5.);
}

TYPED_TEST(MfieldsTest, SetupFields)
{
  using Mfields = TypeParam;

  auto grid = make_grid();
  auto mflds = Mfields{grid, NR_FIELDS, {}};

  setupFields(mflds, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });

  for (int p = 0; p < mflds.n_patches(); ++p) {
    mprintf("p = %d\n", p);

    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
#if 0
	mprintf("[%d, %d, %d] = %06g %06g %06g\n", i, j, k,
		(double) mflds[p](EX, i, j, k), (double) mflds[p](EY, i, j, k), (double) mflds[p](EZ, i, j, k));
#endif
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
