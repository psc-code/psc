
#include "gtest/gtest.h"

#include "testing.hxx"

#include "grid.hxx"
#include "fields.hxx"

#include <limits>

// ---------------------------------------------------------------------------
// DepositNc

template <typename R, typename D>
class DepositNc
{
public:
  using real_t = R;
  using real3_t = Vec3<R>;
  using dim_t = D;

  template <typename F>
  void operator()(F& flds, real3_t x)
  {
    Int3 l;
    real3_t h;
    for (int d = 0; d < 3; d++) {
      l[d] = fint(x[d]);
      h[d] = x[d] - l[d];
    }
    (*this)(flds, l, h);
  }

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h, dim_yz tag)
  {
    flds(0, l[1] + 0, l[2] + 0) += (1. - h[1]) * (1. - h[2]);
    flds(0, l[1] + 0, l[2] + 1) += (1. - h[1]) * h[2];
    flds(0, l[1] + 1, l[2] + 0) += h[1] * (1. - h[2]);
    flds(0, l[1] + 1, l[2] + 1) += h[1] * h[2];
  }

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h)
  {
    (*this)(flds, l, h, dim_t{});
  }
};

// ---------------------------------------------------------------------------

class CalcDivNc
{
public:
  template <typename F>
  auto operator()(const F& flds)
  {
    return (flds.view(_all, _s(1, _), _s(1, _), 1) -
            flds.view(_all, _s(0, -1), _s(1, _), 1)) +
           (flds.view(_all, _s(1, _), _s(1, _), 2) -
            flds.view(_all, _s(1, _), _s(0, -1), 2));
  }
};

// ---------------------------------------------------------------------------

template <typename T>
struct DepositTest : ::testing::Test
{
  using real_t = typename T::real_t;
  using real3_t = Vec3<real_t>;
  using dim_t = typename T::dim_t;
  using Current = typename T::Current;

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

    ldims_ = gt::shape_type<3>(grid().ldims);
  }

  const Grid_t& grid() const { return *grid_; }

  template <typename F>
  void test_charge(real3_t x, const F& rho_ref)
  {
    const real_t eps = std::numeric_limits<real_t>::epsilon();

    DepositNc<real_t, dim_t> deposit;

    auto flds = gt::zeros<real_t>(ldims_);
    deposit(flds, x);

    EXPECT_LT(gt::norm_linf(flds.view(0, _all, _all) - rho_ref), eps)
      << flds.view(0, _all, _all);
  }

  auto calc_current(real3_t xm, real3_t xp, real3_t vxi)
  {
    using fields_view_t = SArrayView<real_t, gt::space::host>;
    using fields_t = curr_cache_t<fields_view_t, dim_xyz>;

    const Grid_t& grid = this->grid();
    Current curr{grid};

    auto flds =
      gt::zeros<real_t>(gt::shape(ldims_[0], ldims_[1], ldims_[2], NR_FIELDS));

    // FIXME, to_kernel() is a hack to get a gtensor_view
    auto flds_view = fields_view_t({0, 0, 0}, flds.to_kernel());
    fields_t J{flds_view};

    real_t qni_wni = 1.;
    Int3 lf = {}, lg = {}; // FIXME, those aren't actually used at all in
                           // Current1VbSplit, at least
    curr.calc_j(J, xm, xp, lf, lg, qni_wni, vxi);
    return flds;
  }

  template <typename F>
  void check_continuity(real3_t xm, real3_t xp, const F& flds)
  {
    const real_t eps = 2. * std::numeric_limits<real_t>::epsilon();

    const Grid_t& grid = this->grid();
    auto rho_m = gt::zeros<real_t>(ldims_);
    auto rho_p = gt::zeros<real_t>(ldims_);
    DepositNc<real_t, dim_t> deposit;
    deposit(rho_m, xm);
    deposit(rho_p, xp);
    auto d_rho = (rho_p - rho_m).view(_all, _s(1, _), _s(1, _));

    CalcDivNc div;
    auto div_j = div(flds.view(_all, _all, _all, _s(JXI, JXI + 3)));
    // std::cout << "d_rho\n" << d_rho.view(0, _all, _all) << "\n";
    // std::cout << "div_j\n" << div_j.view(0, _all, _all) << "\n";
    EXPECT_LT(gt::norm_linf(d_rho + div_j), eps);
  }

  template <typename F>
  void test_current(real3_t xm, real3_t xp, real3_t vxi, const F& jxi_ref,
                    const F& jyi_ref, const F& jzi_ref)
  {
    const real_t eps = 2. * std::numeric_limits<real_t>::epsilon();

    auto flds = calc_current(xm, xp, vxi);

    EXPECT_LT(gt::norm_linf(flds.view(0, _all, _all, JXI) - jxi_ref), eps)
      << flds.view(0, _all, _all, JXI);
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, _all, JYI) - jyi_ref), eps)
      << flds.view(0, _all, _all, JYI);
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, _all, JZI) - jzi_ref), eps)
      << flds.view(0, _all, _all, JZI);

    check_continuity(xm, xp, flds);
  }

  void test_current(real3_t xm, real3_t xp, real3_t vxi)
  {
    auto flds = calc_current(xm, xp, vxi);
    check_continuity(xm, xp, flds);
  }

  Int3 ibn = {0, 0, 0};
  std::unique_ptr<Grid_t> grid_;
  gt::shape_type<3> ldims_;
};

template <typename R, typename D,
          template <typename, typename, typename> class C>
class DepositTestConfig
{
public:
  using real_t = R;
  using dim_t = D;
  using fields_view_t = SArrayView<real_t, gt::space::host>;
  using fields_t = curr_cache_t<fields_view_t, dim_xyz>;
  using Current = C<opt_order_1st, dim_yz, fields_t>;
};

using DepositTestTypes =
  ::testing::Types< // DepositTestConfig<float, dim_yz, Current1vbSplit>,
    DepositTestConfig<double, dim_yz, Current1vbSplit>,
    DepositTestConfig<double, dim_yz, CurrentZigzag>>;

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

TYPED_TEST(DepositTest, CurrentNotMoving)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1., 1.}, xp = {.5, 1., 1.};
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
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, CurrentY)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1., 1.}, xp = {.5, 1.2, 1.};
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
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, CurrentYShift)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1., 1.7}, xp = {.5, 1.2, 1.7};
  real3_t vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0., 0.},
                                    {0., 0.06, 0., 0.},
                                    {0., 0.14, 0., 0.},
                                    {0., 0.  , 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, CurrentZ)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1., 1.3}, xp = {.5, 1., 1.6};
  real3_t vxi = {0., 0., 3.};
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
                                    {0., 0.3, 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, CurrentYCross)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1.9, 1.}, xp = {.5, 2.1, 1.};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  real3_t vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0. , 0., 0.},
                                    {0., 0.1, 0.1, 0.},
                                    {0., 0. , 0., 0.},
                                    {0., 0. , 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, CurrentYCrossShift)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1.9, 1.3}, xp = {.5, 2.1, 1.3};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  real3_t vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0.  , 0.},
                                    {0., 0.07, 0.07, 0.},
                                    {0., 0.03, 0.03, 0.},
                                    {0., 0.  , 0.  , 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, CurrentXY)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1.2, 1.1}, xp = {.5, 1.4, 1.4};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  real3_t vxi = {0., .2, .3};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0., 0.},
                                    {0., 0.15, 0., 0.},
                                    {0., 0.05, 0., 0.},
                                    {0., 0.  , 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0.  , 0.  , 0.},
                                    {0., 0.21, 0.09, 0.},
                                    {0., 0.  , 0.  , 0.},
                                    {0., 0.  , 0.  , 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, CurrentXYCrossShift)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1.9, 1.3}, xp = {.5, 2.1, 1.4};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  real3_t vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0.  , 0.},
                                    {0., 0.0675, 0.0625, 0.},
                                    {0., 0.0325, 0.0375, 0.},
                                    {0., 0.  , 0.  , 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0.0025, 0.095, 0.0025},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(DepositTest, CurrentXYCrossXY)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;

  real3_t xm = {.5, 1.9, 1.6}, xp = {.5, 2.1, 2.2};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  real3_t vxi = {0., .2, 0.};
  if (std::is_same<typename TypeParam::Current,
                   Current1vbSplit<opt_order_1st, dim_yz,
                                   typename TypeParam::fields_t>>::value) {
    // clang-format off
    gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                      {0., 0., 0., 0.},
                                      {0., 0., 0., 0.},
                                      {0., 0., 0., 0.}};
    gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0.  , 0.},
                                      {0., 0.025, 1./600., 0.},
                                      {0., 0.075, .09+1./600, 0.},
                                      {0., 0., 2./300., 0.}};
    gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                      {0., 0.015, 0.38+1/300., 1./600.},
                                      {0., 0., 0.18 + 2./300., 4./300.},
                                      {0., 0., 0., 0.}};
    // clang-format on
    this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
  } else { // CurrentZigzag
    this->test_current(xm, xp, vxi);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
