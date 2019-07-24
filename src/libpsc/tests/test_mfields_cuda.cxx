
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

#ifdef USE_CUDA

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
class TestMfieldsCuda : public ::testing::Test
{};

using TestTypes = ::testing::Types<MfieldsCuda>;

TYPED_TEST_SUITE(TestMfieldsCuda, TestTypes);

TYPED_TEST(TestMfieldsCuda, HostMirror)
{
  using Mfields = TypeParam;

  auto grid = make_grid();
  auto mflds = Mfields{grid, NR_FIELDS, {}};

  setupFields(grid, mflds, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });

  auto h_mflds = hostMirror(mflds);
  copy(mflds, h_mflds);
  
  for (int p = 0; p < mflds.n_patches(); ++p) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	EXPECT_EQ(mflds[p](EX, i, j, k), h_mflds[p](EX, i, j, k));
#if 0
	mprintf("[%d, %d, %d] = %06g %06g\n", i, j, k,
		(double) mflds[p](EX, i, j, k), (double) h_mflds[p](EX, i, j, k));
#endif
    });
  }
}

#endif

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
