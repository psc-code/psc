
#pragma once

// ======================================================================

template <typename Order, typename Dim, typename _fields_t>
struct CurrentZigzag
{
  static_assert(std::is_same<Order, opt_order_1st>::value,
                "CurrentZigzag only works with 1st order");

  using fields_t = _fields_t;
  using real_t = typename fields_t::value_type;
  using Real3 = Vec3<real_t>;

  CurrentZigzag(const Grid_t& grid)
    : dt_(grid.dt),
      dxi_{Real3{1., 1., 1.} / Real3(grid.domain.dx)},
      deposition_(real_t(grid.norm.fnqs / grid.dt) * Real3(grid.domain.dx))
  {}

  void calc_j2_one_cell(fields_t curr_cache, real_t qni_wni, real_t xm[3],
                        real_t xp[3])
  {
    real_t dx[3] = {xp[0] - xm[0], xp[1] - xm[1], xp[2] - xm[2]};
    real_t xa[3] = {.5f * (xm[0] + xp[0]), .5f * (xm[1] + xp[1]),
                    .5f * (xm[2] + xp[2])};

    int i[3];
    for (int d = 0; d < 3; d++) {
      i[d] = fint(xa[d]);
      xa[d] -= i[d];
    }

    deposition_(curr_cache, i, qni_wni, dx, xa, Dim{});
  }

  // ----------------------------------------------------------------------
  // calc_j

  void calc_j(fields_t curr_cache, real_t* xm, real_t* xp, int* lf, int* lg,
              real_t qni_wni, real_t* vxi)
  {
    // this way, we guarantee that the average position in invar directions
    // will remain in the 0th cell
    if (Dim::InvarX::value) {
      xm[0] = .5f;
      xp[0] = xm[0] + vxi[0] * dt_ * dxi_[0];
    }
    if (Dim::InvarY::value) {
      xm[1] = .5f;
      xp[1] = xm[1] + vxi[1] * dt_ * dxi_[1];
    }
    if (Dim::InvarZ::value) {
      xm[1] = .5f;
      xp[1] = xm[1] + vxi[1] * dt_ * dxi_[1];
    }

    real_t xr[3];
    bool crossed = false;
    for (int d = 0; d < 3; d++) {
      int im = fint(xm[d]), ip = fint(xp[d]);
      if (im == ip) {
        xr[d] = .5 * (xm[d] + xp[d]);
      } else {
        xr[d] = std::max(im, ip);
        crossed = true;
      }
    }
    if (crossed) { // could choose this path always
      calc_j2_one_cell(curr_cache, qni_wni, xm, xr);
      calc_j2_one_cell(curr_cache, qni_wni, xr, xp);
    } else {
      calc_j2_one_cell(curr_cache, qni_wni, xm, xp);
    }
  }

private:
  real_t dt_;
  Real3 dxi_;
  psc::CurrentDeposition1vb<fields_t> deposition_;
};
