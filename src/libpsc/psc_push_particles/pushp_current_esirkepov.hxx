
#pragma once

// ----------------------------------------------------------------------
// find_l_minmax

template<typename order>
static void find_l_minmax(int *l1min, int *l1max, int k1, int lg1);

template<>
void find_l_minmax<opt_order_1st>(int *l1min, int *l1max, int k1, int lg1)
{
  if (k1 == lg1) {
    *l1min = 0; *l1max = +1;
  } else if (k1 == lg1 - 1) {
    *l1min = -1; *l1max = +1;
  } else { // (k1 == lg1 + 1)
    *l1min = 0; *l1max = +2;
  }
}

template<>
void find_l_minmax<opt_order_2nd>(int *l1min, int *l1max, int k1, int lg1)
{
  if (k1 == lg1) {
    *l1min = -1; *l1max = +1;
  } else if (k1 == lg1 - 1) {
    *l1min = -2; *l1max = +1;
  } else { // (k1 == lg1 + 1)
    *l1min = -1; *l1max = +2;
  }
}

// ----------------------------------------------------------------------
// Rho1d: charge density

template<typename real_t, typename order>
class Rho1d;

template<typename real_t>
class Rho1d<real_t, opt_order_1st>
{
public:
  static const int N_RHO = 4;
  static const int S_OFF = 1;

  real_t  operator[](int i) const { return s_[i + S_OFF]; }
  real_t& operator[](int i)       { return s_[i + S_OFF]; }

  void zero()
  {
    for (int i = 0; i < N_RHO; i++) {
      s_[i] = 0.;
    }
  }

  void set(int shift, struct ip_coeff_1st<real_t> gg)
{
  (*this)[shift  ] = gg.v0;
  (*this)[shift+1] = gg.v1;
}

private:
  real_t s_[N_RHO];
};

template<typename real_t>
class Rho1d<real_t, opt_order_2nd>
{
public:
  static const int N_RHO = 5;
  static const int S_OFF = 2;

  real_t  operator[](int i) const { return s_[i + S_OFF]; }
  real_t& operator[](int i)       { return s_[i + S_OFF]; }

  void zero()
  {
    for (int i = 0; i < N_RHO; i++) {
      s_[i] = 0.;
    }
  }

  void set(int shift, struct ip_coeff_2nd<real_t> gg)
{
  // FIXME: It appears that gm/g0/g1 can be used instead of what's calculated here
  // but it needs checking.
  real_t h = gg.h;
  (*this)[shift-1] = .5f * (1.5f-std::abs(h-1.f)) * (1.5f-std::abs(h-1.f));
  (*this)[shift  ] = .75f - std::abs(h) * std::abs(h);
  (*this)[shift+1] = .5f * (1.5f-std::abs(h+1.f)) * (1.5f-std::abs(h+1.f));
}

private:
  real_t s_[N_RHO];
};

template<typename order, typename IP, typename Invar>
struct CurrentDir
{
  using real_t = typename IP::real_t;
  using Rho1d_t = Rho1d<real_t, order>;
  using ip_coeff_t = typename IP::ip_coeff_t;

  void charge_before(ip_coeff_t g)
  {
    lg = g.l;
    s0.set(0, g);
  };

  void charge_after(real_t xx)
  {
    ip_coeff_t g;
    g.set(xx);
    k = g.l;
    s1.zero();
    s1.set(k - lg, g);
  }

  void prep(real_t qni_wni, real_t vv, real_t fnqxyzs, real_t fnqs)
  {
    for (int i = -s1.S_OFF + 1; i <= 1; i++) {
      s1[i] -= s0[i];
    }
    find_l_minmax<order>(&lmin, &lmax, k, lg);
    fnq = qni_wni * fnqxyzs;
  }

  Rho1d_t s0 = {}, s1;
  int lg, k;
  int lmin, lmax;
  real_t fnq;
};

// true_type -> this is an invariant direction
template<typename order, typename IP>
struct CurrentDir<order, IP, std::true_type>
{
  using real_t = typename IP::real_t;
  using ip_coeff_t = typename IP::ip_coeff_t;

  void charge_before(ip_coeff_t g) {}

  void charge_after(real_t xx) {}

  void prep(real_t qni_wni, real_t vv, real_t fnqxyzs, real_t fnqs)
  {
    fnqv = vv * qni_wni * fnqs;
  }

  real_t fnqv;
};

// ======================================================================
// Current

template<typename order, typename dim, typename IP, typename Fields>
struct Current
{
  using real_t = typename IP::real_t;
  using Real3 = Vec3<real_t>;

  Current(const Grid_t& grid)
    : dxi_{ Real3{1., 1. , 1.} / Real3{grid.domain.dx} },
      fnqs_(grid.norm.fnqs)
  {
    fnqxs_ = grid.domain.dx[0] * fnqs_ / grid.dt;
    fnqys_ = grid.domain.dx[1] * fnqs_ / grid.dt;
    fnqzs_ = grid.domain.dx[2] * fnqs_ / grid.dt;
  }

  void charge_before(const IP& ip)
  {
    x.charge_before(ip.cx.g);
    y.charge_before(ip.cy.g);
    z.charge_before(ip.cz.g);
  }

  void charge_after(real_t xx[3])
  {
    x.charge_after(xx[0] * dxi_[0]);
    y.charge_after(xx[1] * dxi_[1]);
    z.charge_after(xx[2] * dxi_[2]);
  }

  void prep(real_t qni_wni, real_t vv[3])
  {
    x.prep(qni_wni, vv[0], fnqxs_, fnqs_);
    y.prep(qni_wni, vv[1], fnqys_, fnqs_);
    z.prep(qni_wni, vv[2], fnqzs_, fnqs_);
  }

  void calc(Fields& J)
  {
    // Dispatch by function overloading, which needs to pass order and dim as tags
    calc(order{}, dim{}, J);
  }

  void calc(opt_order_2nd o, dim_1 d, Fields& J);
  void calc(opt_order_2nd o, dim_y d, Fields& J);
  void calc(opt_order_2nd o, dim_z d, Fields& J);
  void calc(opt_order_2nd o, dim_xy d, Fields& J);
  void calc(opt_order_2nd o, dim_xz d, Fields& J);
  void calc(opt_order_1st o, dim_xz d, Fields& J);
  void calc(opt_order_1st o, dim_yz d, Fields& J);
  void calc(opt_order_2nd o, dim_yz d, Fields& J);
  void calc(opt_order_2nd o, dim_xyz d, Fields& J);

  CurrentDir<order, IP, typename dim::InvarX> x;
  CurrentDir<order, IP, typename dim::InvarY> y;
  CurrentDir<order, IP, typename dim::InvarZ> z;

private:
  Real3 dxi_;
  real_t fnqs_;
  real_t fnqxs_, fnqys_, fnqzs_;
};

// ======================================================================

template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_2nd o, dim_1 d, Fields& J)
{
  real_t jxh = x.fnqv;
  real_t jyh = y.fnqv;
  real_t jzh = z.fnqv;

  J(JXI, 0,0,0) += jxh;
  J(JYI, 0,0,0) += jyh;
  J(JZI, 0,0,0) += jzh;
};

template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_2nd o, dim_y d, Fields& J)
{
  real_t jyh = 0.f;

  for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
    real_t wx = y.s0[l2] + .5f * y.s1[l2];
    real_t wy = y.s1[l2];
    real_t wz = y.s0[l2] + .5f * y.s1[l2];

    real_t jxh = x.fnqv*wx;
    jyh -= y.fnq*wy;
    real_t jzh = z.fnqv*wz;

    J(JXI, 0,y.lg+l2,0) += jxh;
    J(JYI, 0,y.lg+l2,0) += jyh;
    J(JZI, 0,y.lg+l2,0) += jzh;
  }
}

template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_2nd o, dim_z d, Fields& J)
{
  real_t jzh = 0.f;
  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    real_t wx = z.s0[l3] + .5f * z.s1[l3];
    real_t wy = z.s0[l3] + .5f * z.s1[l3];
    real_t wz = z.s1[l3];

    real_t jxh = x.fnqv*wx;
    real_t jyh = y.fnqv*wy;
    jzh -= z.fnq*wz;

    J(JXI, 0,0,z.lg+l3) += jxh;
    J(JYI, 0,0,z.lg+l3) += jyh;
    J(JZI, 0,0,z.lg+l3) += jzh;
  }
};

template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_2nd o, dim_xy d, Fields& J)
{
  for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
    real_t jxh = 0.f;
    for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
      real_t wx = x.s1[l1] * (y.s0[l2] + .5f*y.s1[l2]);
      real_t wz = x.s0[l1] * y.s0[l2]
	+ .5f * x.s1[l1] * y.s0[l2]
	+ .5f * x.s0[l1] * y.s1[l2]
	+ (1.f/3.f) * x.s1[l1] * y.s1[l2];

      jxh -= x.fnq*wx;
      J(JXI, x.lg+l1,y.lg+l2,0) += jxh;
      J(JZI, x.lg+l1,y.lg+l2,0) += z.fnqv * wz;
    }
  }
  for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
    real_t jyh = 0.f;
    for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
      real_t wy = y.s1[l2] * (x.s0[l1] + .5f*x.s1[l1]);

      jyh -= y.fnq*wy;
      J(JYI, x.lg+l1,y.lg+l2,0) += jyh;
    }
  }
}

// FIXME, 1st/2d duplicated
template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_1st o, dim_xz d, Fields& J)
{
  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    real_t jxh = 0.f;
    for (int l1 = x.lmin; l1 < x.lmax; l1++) {
      real_t wx = x.s1[l1] * (z.s0[l3] + .5f*z.s1[l3]);
      jxh -= x.fnq*wx;
      J(JXI, x.lg+l1,0,z.lg+l3) += jxh;
    }
  }

  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
      real_t wy = x.s0[l1] * z.s0[l3]
	+ .5f * x.s1[l1] * z.s0[l3]
	+ .5f * x.s0[l1] * z.s1[l3]
	+ (1.f/3.f) * x.s1[l1] * z.s1[l3];
      real_t jyh = y.fnqv * wy;
      J(JYI, x.lg+l1,0,z.lg+l3) += jyh;
    }
  }
  for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
    real_t jzh = 0.f;
    for (int l3 = z.lmin; l3 < z.lmax; l3++) {
      real_t wz = z.s1[l3] * (x.s0[l1] + .5f*x.s1[l1]);
      jzh -= z.fnq*wz;
      J(JZI, x.lg+l1,0,z.lg+l3) += jzh;
    }
  }
}

template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_2nd o, dim_xz d, Fields& J)
{
  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    real_t jxh = 0.f;
    for (int l1 = x.lmin; l1 < x.lmax; l1++) {
      real_t wx = x.s1[l1] * (z.s0[l3] + .5f*z.s1[l3]);
      jxh -= x.fnq*wx;
      J(JXI, x.lg+l1,0,z.lg+l3) += jxh;
    }
  }

  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
      real_t wy = x.s0[l1] * z.s0[l3]
	+ .5f * x.s1[l1] * z.s0[l3]
	+ .5f * x.s0[l1] * z.s1[l3]
	+ (1.f/3.f) * x.s1[l1] * z.s1[l3];
      real_t jyh = y.fnqv * wy;
      J(JYI, x.lg+l1,0,z.lg+l3) += jyh;
    }
  }
  for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
    real_t jzh = 0.f;
    for (int l3 = z.lmin; l3 < z.lmax; l3++) {
      real_t wz = z.s1[l3] * (x.s0[l1] + .5f*x.s1[l1]);
      jzh -= z.fnq*wz;
      J(JZI, x.lg+l1,0,z.lg+l3) += jzh;
    }
  }
}

template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_1st o, dim_yz d, Fields& J)
{
  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
      real_t wx = y.s0[l2] * z.s0[l3]
	+ .5f * y.s1[l2] * z.s0[l3]
	+ .5f * y.s0[l2] * z.s1[l3]
	+ (1.f/3.f) * y.s1[l2] * z.s1[l3];
      real_t jxh = x.fnqv * wx;
      J(JXI, 0,y.lg+l2,z.lg+l3) += jxh;
    }
  }

  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    real_t jyh = 0.f;
    for (int l2 = y.lmin; l2 < y.lmax; l2++) {
      real_t wy = y.s1[l2] * (z.s0[l3] + .5f*z.s1[l3]);
      jyh -= y.fnq*wy;
      J(JYI, 0,y.lg+l2,z.lg+l3) += jyh;
    }
  }

  for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
    real_t jzh = 0.f;
    for (int l3 = z.lmin; l3 < z.lmax; l3++) {
      real_t wz = z.s1[l3] * (y.s0[l2] + .5f*y.s1[l2]);
      jzh -= z.fnq*wz;
      J(JZI, 0,y.lg+l2,z.lg+l3) += jzh;
    }
  }
}

#define JZH(i) jzh[i+2]
template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_2nd o, dim_yz d, Fields& J)
{
  real_t jxh;
  real_t jyh;
  real_t jzh[5];

  for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
    JZH(l2) = 0.f;
  }
  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    jyh = 0.f;
    for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
      real_t wx = y.s0[l2] * z.s0[l3]
	+ .5f * y.s1[l2] * z.s0[l3]
	+ .5f * y.s0[l2] * z.s1[l3]
	+ (1.f/3.f) * y.s1[l2] * z.s1[l3];
      real_t wy = y.s1[l2] * (z.s0[l3] + .5f*z.s1[l3]);
      real_t wz = z.s1[l3] * (y.s0[l2] + .5f*y.s1[l2]);

      jxh = x.fnqv*wx;
      jyh -= y.fnq*wy;
      JZH(l2) -= z.fnq*wz;

      J(JXI, 0,y.lg+l2,z.lg+l3) += jxh;
      J(JYI, 0,y.lg+l2,z.lg+l3) += jyh;
      J(JZI, 0,y.lg+l2,z.lg+l3) += JZH(l2);
    }
  }
}

template<typename order, typename dim, typename IP, typename Fields>
void Current<order, dim, IP, Fields>::calc(opt_order_2nd o, dim_xyz d, Fields& J)
{
  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
      real_t jxh = 0.f;
      for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
	real_t wx = x.s1[l1] * (y.s0[l2] * z.s0[l3] +
					   .5f * y.s1[l2] * z.s0[l3] +
					   .5f * y.s0[l2] * z.s1[l3] +
					   (1.f/3.f) * y.s1[l2] * z.s1[l3]);

	jxh -= x.fnq*wx;
	J(JXI, x.lg+l1,y.lg+l2,z.lg+l3) += jxh;
      }
    }
  }

  for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
    for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
      real_t jyh = 0.f;
      for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
	real_t wy = y.s1[l2] * (x.s0[l1] * z.s0[l3] +
					   .5f * x.s1[l1] * z.s0[l3] +
					   .5f * x.s0[l1] * z.s1[l3] +
					   (1.f/3.f) * x.s1[l1]*z.s1[l3]);

	jyh -= y.fnq*wy;
	J(JYI, x.lg+l1,y.lg+l2,z.lg+l3) += jyh;
      }
    }
  }

  for (int l2 = y.lmin; l2 <= y.lmax; l2++) {
    for (int l1 = x.lmin; l1 <= x.lmax; l1++) {
      real_t jzh = 0.f;
      for (int l3 = z.lmin; l3 <= z.lmax; l3++) {
	real_t wz = z.s1[l3] * (x.s0[l1] * y.s0[l2] +
					   .5f * x.s1[l1] * y.s0[l2] +
					   .5f * x.s0[l1] * y.s1[l2] +
					   (1.f/3.f) * x.s1[l1]*y.s1[l2]);

	jzh -= z.fnq*wz;
	J(JZI, x.lg+l1,y.lg+l2,z.lg+l3) += jzh;
      }
    }
  }
 }

