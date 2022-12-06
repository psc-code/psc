
#include "gtest/gtest.h"

#include "testing.hxx"

#include "grid.hxx"
#include "fields.hxx"
#include "psc_particles_double.h"
#include "psc_fields_c.h"
#include "push_particles.hxx"
#include "setup_fields.hxx"
#include "setup_particles.hxx"

// ---------------------------------------------------------------------------
// DepositNc

template <typename R>
class DepositNc
{
public:
  using real_t = R;
  using real3_t = Vec3<R>;

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h)
  {
    flds(0, l[1] + 0, l[2] + 0) += (1. - h[1]) * (1. - h[2]);
    flds(0, l[1] + 0, l[2] + 1) += (1. - h[1]) * h[2];
    flds(0, l[1] + 1, l[2] + 0) += h[1] * (1. - h[2]);
    flds(0, l[1] + 1, l[2] + 1) += h[1] * h[2];
  }
};

// ---------------------------------------------------------------------------

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
  void test_charge(real3_t x, const F& rho_ref)
  {
    const real_t eps = 1e-7;

    Int3 l;
    real3_t h;
    for (int d = 0; d < 3; d++) {
      l[d] = fint(x[d]);
      h[d] = x[d] - l[d];
    }
    DepositNc<real_t> deposit;

    auto flds = gt::zeros<real_t>(gt::shape(1, 4, 4));
    deposit(flds, l, h);

    // std::cout << flds.view(0, _all, _all);
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, _all) - rho_ref), eps);
  }

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

TYPED_TEST(DepositTest, ChargeCenter)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t x = {.5, 1.5, 1.5};
  // clang-format off
  gt::gtensor<real_t, 2> rho_ref = {{0., 0.  , 0.  , 0.},
                                    {0., 0.25, 0.25, 0.},
                                    {0., 0.25, 0.25, 0.},
                                    {0., 0.  , 0.  , 0.}};
  // clang-format on
  this->test_charge(x, rho_ref);
}

TYPED_TEST(DepositTest, ChargeLowerLeft)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t x = {.5, 1., 1.};
  // clang-format off
  gt::gtensor<real_t, 2> rho_ref = {{0., 0., 0., 0.},
                                    {0., 1., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_charge(x, rho_ref);
}

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
