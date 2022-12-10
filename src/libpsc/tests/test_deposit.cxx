
#include "gtest/gtest.h"

#include <psc/deposit.hxx>

#include <gtensor/reductions.h>
#include <mpi.h>

#include <limits>

using namespace gt::placeholders;

// ---------------------------------------------------------------------------

template <typename T>
struct DepositTest : ::testing::Test
{
  using real_t = typename T::real_t;
  using real3_t = gt::sarray<real_t, 3>;
  using dim_t = typename T::dim_t;

  DepositTest() : ldims_{4, 4, 4}
  {
    if (dim_t::InvarX::value) {
      ldims_[0] = 1;
    }
  }

  gt::shape_type<3> ldims_;
  const real_t eps = std::numeric_limits<real_t>::epsilon();
};

template <typename R, typename D>
class DepositTestConfig
{
public:
  using real_t = R;
  using dim_t = D;
};

using DepositTestTypes = ::testing::Types<DepositTestConfig<double, dim_yz>,
                                          DepositTestConfig<double, dim_xyz>>;

TYPED_TEST_SUITE(DepositTest, DepositTestTypes);

TYPED_TEST(DepositTest, ChargeCenter)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;
  using dim_t = typename self_type::dim_t;

  real3_t x = {1.5, 1.5, 1.5};
  real_t val = .1;
  auto ibn = gt::shape(0, 0, 0);
  auto rho_ref = gt::zeros<real_t>(this->ldims_ + 2 * ibn);
  if (std::is_same<dim_t, dim_xyz>::value) {
    rho_ref.view(_s(1, 3), _s(1, 3), _s(1, 3)) = val / 8.f;
  } else if (std::is_same<dim_t, dim_yz>::value) {
    rho_ref.view(0, _s(1, 3), _s(1, 3)) = val / 4.f;
  }

  auto flds = gt::zeros_like(rho_ref);
  psc::deposit::nc<dim_t>(flds, -ibn, x, val);
  EXPECT_LT(gt::norm_linf(flds - rho_ref), this->eps) << flds << "\nref\n"
                                                      << rho_ref;
}

TYPED_TEST(DepositTest, ChargeLowerLeft)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;
  using dim_t = typename self_type::dim_t;

  real3_t x = {1., 1., 1.};
  real_t val = .1;
  auto ibn = gt::shape(0, 0, 0);
  auto rho_ref = gt::zeros<real_t>(this->ldims_ + 2 * ibn);
  if (std::is_same<dim_t, dim_xyz>::value) {
    rho_ref(1, 1, 1) = val;
  } else if (std::is_same<dim_t, dim_yz>::value) {
    rho_ref(0, 1, 1) = val;
  }

  auto flds = gt::zeros_like(rho_ref);
  psc::deposit::nc<dim_t>(flds, -ibn, x, val);
  EXPECT_LT(gt::norm_linf(flds - rho_ref), this->eps) << flds << "\nref\n"
                                                      << rho_ref;
}

TYPED_TEST(DepositTest, ChargeCenterWithBnd)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;
  using dim_t = typename self_type::dim_t;

  real3_t x = {.5, .5, .5};
  real_t val = 1.;
  auto ibn = gt::shape(1, 1, 1);
  auto rho_ref = gt::zeros<real_t>(this->ldims_ + 2 * ibn);
  if (std::is_same<dim_t, dim_xyz>::value) {
    rho_ref.view(_s(1, 3), _s(1, 3), _s(1, 3)) = val / 8.f;
  } else if (std::is_same<dim_t, dim_yz>::value) {
    rho_ref.view(0, _s(1, 3), _s(1, 3)) = val / 4.f;
  }

  auto flds = gt::zeros_like(rho_ref);
  psc::deposit::nc<dim_t>(flds, -ibn, x, val);
  EXPECT_LT(gt::norm_linf(flds - rho_ref), this->eps) << flds << "\nref\n"
                                                      << rho_ref;
}

TYPED_TEST(DepositTest, ChargeCenterWithBndCc)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;
  using dim_t = typename self_type::dim_t;

  real3_t x = {.5, .5, .5};
  real_t val = 1.;
  auto ibn = gt::shape(1, 1, 1);
  auto rho_ref = gt::zeros<real_t>(this->ldims_ + 2 * ibn);
  if (std::is_same<dim_t, dim_xyz>::value) {
    rho_ref(1, 1, 1) = val;
  } else if (std::is_same<dim_t, dim_yz>::value) {
    rho_ref(0, 1, 1) = val;
  }

  auto flds = gt::zeros_like(rho_ref);
  psc::deposit::cc<dim_t>(flds, -ibn, x, val);
  EXPECT_LT(gt::norm_linf(flds - rho_ref), this->eps) << flds << "\nref\n"
                                                      << rho_ref;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
