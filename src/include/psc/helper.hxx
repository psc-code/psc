
#pragma once

namespace psc
{
namespace helper
{

template <typename E1, typename E2>
void print_diff(const E1& e1, const E2& e2, double eps)
{
  auto&& h_e1 = gt::host_mirror(e1);
  auto&& h_e2 = gt::host_mirror(e2);
  gt::copy(gt::eval(e1), h_e1);
  gt::copy(gt::eval(e2), h_e2);
  gt::launch_host<5>(h_e1.shape(), [=](int i, int j, int k, int m, int p) {
    auto val_e1 = h_e1(i, j, k, m, p);
    auto val_e2 = h_e2(i, j, k, m, p);
    if (std::abs(val_e1 + val_e2) > eps) {
      mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", p, i, j, k, val_e1, val_e2,
              val_e1 - val_e2);
    }
  });
}

template <typename E1, typename E2>
void print_diff_3d(const E1& e1, const E2& e2, double eps)
{
  auto&& h_e1 = gt::host_mirror(e1);
  auto&& h_e2 = gt::host_mirror(e2);
  gt::copy(gt::eval(e1), h_e1);
  gt::copy(gt::eval(e2), h_e2);
  gt::launch_host<3>(h_e1.shape(), [=](int i, int j, int k) {
    auto val_e1 = h_e1(i, j, k);
    auto val_e2 = h_e2(i, j, k);
    if (std::abs(val_e1 + val_e2) > eps) {
      mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", i, j, k, val_e1, val_e2,
              val_e1 - val_e2);
    }
  });
}

} // namespace helper
} // namespace psc
