
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
  auto h_mflds_ref = MfieldsSingle{grid, NR_FIELDS, {}};

  setupFields(mflds, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });
  setupFields(h_mflds_ref, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });

  auto&& gt = mflds.storage();
  auto&& h_gt = gt::host_mirror(gt);
  gt::copy(gt, h_gt);
  auto&& h_gt_ref = h_mflds_ref.storage();

  gt::launch<5, gt::space::host>(
    h_gt.shape(), [=](int i, int j, int k, int m, int p) {
      EXPECT_EQ(h_gt(i, j, k, m, p), h_gt_ref(i, j, k, m, p));
    });
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
