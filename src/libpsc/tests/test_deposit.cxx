
#include "gtest/gtest.h"

#include "testing.hxx"

#include "grid.hxx"
#include "fields.hxx"
#include "psc_particles_double.h"
#include "psc_fields_c.h"
#include "push_particles.hxx"
#include "setup_fields.hxx"
#include "setup_particles.hxx"

template <typename T>
struct DepositTest : ::testing::Test
{
  using real_t = typename T::Mfields::real_t;
  using real3_t = Vec3<real_t>;

  DepositTest()
  {
    double dt = 1.;
    Int3 gdims = {1, 4, 4};
    Double3 length = {1., 4., 4.};
    psc::grid::Domain<double> grid_domain{gdims, length};
    psc::grid::BC grid_bc;
    Grid_t::Normalization grid_norm;
    std::vector<Grid_t::Kind> kinds;

    grid_.reset(
      new Grid_t{grid_domain, grid_bc, kinds, grid_norm, dt, 1, this->ibn});
  }

  const Grid_t& grid() const { return *grid_; }

  template <typename F>
  void test_deposit(real3_t xm, real3_t xp, Int3 lf, Int3 lg, real3_t vxi,
                    const F& jxi_ref, const F& jyi_ref, const F& jzi_ref)
  {
    using Mfields = typename T::Mfields;
    using fields_t = curr_cache_t<typename Mfields::fields_view_t, dim_xyz>;
    using Current = Current1vbSplit<opt_order_1st, dim_yz, fields_t>;

    const real_t eps = 1e-7;

    Current curr{this->grid()};

    Mfields mflds{this->grid(), NR_FIELDS, this->ibn};
    int p = 0;
    auto flds = mflds[p];
    ASSERT_EQ(gt::norm_linf(flds.storage()), 0.);

    typename Current::fields_t J{flds};

    real_t qni_wni = 1.;
    curr.calc_j(J, xm, xp, lf, lg, qni_wni, vxi);

    // std::cout << "JXI\n" << flds.storage().view(0, _all, _all, JXI) << "\n";
    // std::cout << "JYI\n" << flds.storage().view(0, _all, _all, JYI) << "\n";
    // std::cout << "JZI\n" << flds.storage().view(0, _all, _all, JZI) << "\n";

    EXPECT_LT(gt::norm_linf(flds.storage().view(0, _all, _all, JXI) - jxi_ref),
              eps);
    EXPECT_LT(gt::norm_linf(flds.storage().view(0, _all, _all, JYI) - jyi_ref),
              eps);
    EXPECT_LT(gt::norm_linf(flds.storage().view(0, _all, _all, JZI) - jzi_ref),
              eps);
  }

  Int3 ibn = {0, 0, 0};

  std::unique_ptr<Grid_t> grid_;
};

using DepositTestTypes = ::testing::Types<TestConfig1vbec3dSingleYZ>;

TYPED_TEST_SUITE(DepositTest, DepositTestTypes);

TYPED_TEST(DepositTest, NotMoving)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1., 1.}, xp = {.5, 1., 1.};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  real3_t vxi = {0., 0., 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_deposit(xm, xp, lf, lg, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, MovingY)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1., 1.}, xp = {.5, 1.2, 1.};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  real3_t vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0. , 0., 0.},
                                    {0., 0.2, 0., 0.},
                                    {0., 0. , 0., 0.},
                                    {0., 0. , 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_deposit(xm, xp, lf, lg, vxi, jxi_ref, jyi_ref, jzi_ref);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
