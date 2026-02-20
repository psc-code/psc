
#pragma once

namespace psc
{
namespace helper
{

namespace detail
{

template <typename E1, typename E2, typename Enable = void>
struct print_diff
{
  static void run(const E1& e1, const E2& e2, double eps)
  {
    static_assert(gt::expr_dimension<E1>() == 5, "e1");
    static_assert(gt::expr_dimension<E2>() == 5, "e2");
  };
};

template <typename E1, typename E2>
struct print_diff<E1, E2,
                  std::enable_if_t<gt::expr_dimension<E1>() == 5 &&
                                   gt::expr_dimension<E2>() == 5>>
{
  static void run(const E1& e1, const E2& e2, double eps)
  {
    auto&& d_e1 = gt::eval(e1);
    auto&& d_e2 = gt::eval(e2);
    auto&& h_e1 = gt::host_mirror(d_e1);
    auto&& h_e2 = gt::host_mirror(d_e2);
    gt::copy(d_e1, h_e1);
    gt::copy(d_e2, h_e2);
    gt::launch_host<5>(h_e1.shape(), [=](int i, int j, int k, int m, int p) {
      auto val_e1 = h_e1(i, j, k, m, p);
      auto val_e2 = h_e2(i, j, k, m, p);
      if (std::abs(val_e1 - val_e2) > eps) {
        mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", p, i, j, k, val_e1,
                val_e2, val_e1 - val_e2);
      }
    });
  }
};

template <typename E1, typename E2>
struct print_diff<E1, E2,
                  std::enable_if_t<gt::expr_dimension<E1>() == 3 &&
                                   gt::expr_dimension<E2>() == 3>>
{
  static void run(const E1& e1, const E2& e2, double eps)
  {
    auto&& d_e1 = gt::eval(e1);
    auto&& d_e2 = gt::eval(e2);
    auto&& h_e1 = gt::host_mirror(d_e1);
    auto&& h_e2 = gt::host_mirror(d_e2);
    gt::copy(d_e1, h_e1);
    gt::copy(d_e2, h_e2);
    gt::launch_host<3>(h_e1.shape(), [=](int i, int j, int k) {
      auto val_e1 = h_e1(i, j, k);
      auto val_e2 = h_e2(i, j, k);
      if (std::abs(val_e1 - val_e2) > eps) {
        mprintf("(%d,%d,%d): %g -- %g diff %g\n", i, j, k, val_e1, val_e2,
                val_e1 - val_e2);
      }
    });
  }
};

} // namespace detail

template <typename E1, typename E2>
void print_diff(const E1& e1, const E2& e2, double eps)
{
  detail::print_diff<E1, E2>::run(e1, e2, eps);
}

} // namespace helper
} // namespace psc
