
#pragma once

// J. Villasenor and O. Buneman, "Rigorous charge conservation for local
// electromagnetic field solvers", Computer Physics Communications 69 (1992) 306

// ======================================================================

template <typename Order, typename Dim, typename _fields_t>
struct Current1vbSplit
{
  static_assert(std::is_same<Order, opt_order_1st>::value,
                "Current1vbSplit only works with 1st order");

  using fields_t = _fields_t;
  using real_t = typename fields_t::value_type;
  using Real3 = Vec3<real_t>;

  Current1vbSplit(const Grid_t& grid)
    : dt_(grid.dt),
      dxi_{grid.domain.dx.inv()},
      deposition_(real_t(grid.norm.fnqs / grid.dt) * Real3(grid.domain.dx))
  {}

  void calc_j2_one_cell(fields_t curr_cache, real_t qni_wni, const Real3& xm,
                        const Real3& xp)
  {
    Real3 dx = xp - xm;
    Real3 xa = real_t(0.5) * (xp + xm);
    Int3 i = xa.fint();
    xa -= Real3(i);
    deposition_(curr_cache, i, qni_wni, dx, xa, Dim{});
  }

  static Real3 calc_split_x1(int dim, int im, int ip, const Real3& xm,
                             const Real3& xp)
  {
    // boundary is at the lower edge of the cell with the higher index
    real_t bnd = std::max(im, ip);
    real_t frac = (bnd - xm[dim]) / (xp[dim] - xm[dim]);

    Real3 x1; // = xm + frac * (xp - xm), but want exact x1[d]
    for (int d = 0; d < 3; d++) {
      if (d == dim) {
        x1[d] = bnd;
      } else {
        x1[d] = xm[d] + frac * (xp[d] - xm[d]);
      }
    }
    return x1;
  }

  // TODO c++17 combine these calc_j2_split_dim_* with if constexpr

  void calc_j2_split_dim_x(fields_t curr_cache, real_t qni_wni, const Real3& xm,
                           const Real3& xp)
  {
    constexpr int dim = 0;
    int im = fint(xm[dim]);
    int ip = fint(xp[dim]);
    if (Dim::is_invar(dim) || ip == im) {
      calc_j2_one_cell(curr_cache, qni_wni, xm, xp);
    } else {
      Real3 x1 = calc_split_x1(dim, im, ip, xm, xp);
      calc_j2_one_cell(curr_cache, qni_wni, xm, x1);
      calc_j2_one_cell(curr_cache, qni_wni, x1, xp);
    }
  }

  void calc_j2_split_dim_y(fields_t curr_cache, real_t qni_wni, const Real3& xm,
                           const Real3& xp)
  {
    constexpr int dim = 1;
    int im = fint(xm[dim]);
    int ip = fint(xp[dim]);
    if (Dim::is_invar(dim) || ip == im) {
      calc_j2_split_dim_x(curr_cache, qni_wni, xm, xp);
    } else {
      Real3 x1 = calc_split_x1(dim, im, ip, xm, xp);
      calc_j2_split_dim_x(curr_cache, qni_wni, xm, x1);
      calc_j2_split_dim_x(curr_cache, qni_wni, x1, xp);
    }
  }

  void calc_j2_split_dim_z(fields_t curr_cache, real_t qni_wni, const Real3& xm,
                           const Real3& xp)
  {
    constexpr int dim = 2;
    int im = fint(xm[dim]);
    int ip = fint(xp[dim]);
    if (Dim::is_invar(dim) || ip == im) {
      calc_j2_split_dim_y(curr_cache, qni_wni, xm, xp);
    } else {
      Real3 x1 = calc_split_x1(dim, im, ip, xm, xp);
      calc_j2_split_dim_y(curr_cache, qni_wni, xm, x1);
      calc_j2_split_dim_y(curr_cache, qni_wni, x1, xp);
    }
  }

  // ----------------------------------------------------------------------
  // calc_j

  void calc_j(fields_t curr_cache, Real3& xm, Real3& xp, const Int3& lf,
              const Int3& lg, real_t qni_wni, Real3& vxi)
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
      xm[2] = .5f;
      xp[2] = xm[2] + vxi[2] * dt_ * dxi_[2];
    }

    calc_j2_split_dim_z(curr_cache, qni_wni, xm, xp);
  }

private:
  real_t dt_;
  Real3 dxi_;
  psc::CurrentDeposition1vb<fields_t> deposition_;
};
