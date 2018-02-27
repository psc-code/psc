
#include <mrc_profile.h>

#include "inc_defs.h"

#define CACHE_EM_J 1

using real_t = mparticles_t::real_t;

#include "fields.hxx"
#include "inc_params.c"
#include "inc_push.c"
#include "inc_cache.c"
#include "interpolate.hxx"

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
// charge density

template<typename order>
class Rho1d;

template<>
class Rho1d<opt_order_1st>
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

template<>
class Rho1d<opt_order_2nd>
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
  using Rho1d_t = Rho1d<order>;
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
  using Real3 = Vec3<real_t>;

  Current(const Grid_t& grid)
    : dxi_{ Real3{1., 1. , 1.} / grid.dx }
  {}

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
    x.prep(qni_wni, vv[0], c_prm.fnqxs, c_prm.fnqs);
    y.prep(qni_wni, vv[1], c_prm.fnqys, c_prm.fnqs);
    z.prep(qni_wni, vv[2], c_prm.fnqzs, c_prm.fnqs);
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

template<typename fields_t, typename dim>
struct CacheFields;

template<typename fields_t>
struct CacheFields<fields_t, dim_yz>
{
  static fields_t from_em(fields_t flds)
  {
    fields_t fld = fields_t(flds.ib, flds.im, 9);
    Fields3d<fields_t> F(flds), F_CACHE(fld);
    // FIXME, can do -1 .. 2? NO!, except maybe for 1st order
    // Has to be at least -2 .. +3 because of staggering
    // FIXME, get rid of caching since it's no different from the actual
    // fields...
    for (int iz = fld.ib[2]; iz < fld.ib[2] + fld.im[2]; iz++) {
      for (int iy = fld.ib[1]; iy < fld.ib[1] + fld.im[1]; iy++) {
	F_CACHE(EX, 0,iy,iz) = F(EX, 0,iy,iz);
	F_CACHE(EY, 0,iy,iz) = F(EY, 0,iy,iz);
	F_CACHE(EZ, 0,iy,iz) = F(EZ, 0,iy,iz);
	F_CACHE(HX, 0,iy,iz) = F(HX, 0,iy,iz);
	F_CACHE(HY, 0,iy,iz) = F(HY, 0,iy,iz);
	F_CACHE(HZ, 0,iy,iz) = F(HZ, 0,iy,iz);
      }
    }
    return fld;
  }

  static void to_j(fields_t fld, fields_t flds)
  {
    Fields3d<fields_t> F(flds), F_CACHE(fld);
    for (int iz = fld.ib[2]; iz < fld.ib[2] + fld.im[2]; iz++) {
      for (int iy = fld.ib[1]; iy < fld.ib[1] + fld.im[1]; iy++) {
	F(JXI, 0,iy,iz) += F_CACHE(JXI, 0,iy,iz);
	F(JYI, 0,iy,iz) += F_CACHE(JYI, 0,iy,iz);
	F(JZI, 0,iy,iz) += F_CACHE(JZI, 0,iy,iz);
      }
    }
    fld.dtor();
  }
};

template<typename fields_t, typename dim>
struct CacheFieldsNone
{
  static fields_t from_em(fields_t flds)
  {
    return flds;
  }

  static void to_j(fields_t flds_cache, fields_t flds)
  {}
};

template<typename C>
struct PushParticles__
{
  using Mparticles = typename C::Mparticles;
  using Mfields = typename C::Mfields;
  using fields_t = typename Mfields::fields_t;
  using CacheFields_t = typename C::CacheFields;

  static void push_mprts(Mparticles& mprts, Mfields& mflds)
  {
    static int pr;
    if (!pr) {
      pr = prof_register(__func__, 1., 0, 0);
    }

    prof_start(pr);
    for (int p = 0; p < mprts.n_patches(); p++) {
      // FIXME, in the cache case can't we just skip this and just set j when copying back?
      mflds[p].zero(JXI, JXI + 3);
      CacheFields_t cache;
      auto flds = cache.from_em(mflds[p]);
      do_push_part(flds, mprts[p]);
      cache.to_j(flds, mflds[p]);
    }
    prof_stop(pr);
  }

private:
  static void do_push_part(fields_t flds, typename Mparticles::patch_t& prts)
  {
    using dim = typename C::dim;
    using IP = InterpolateEM<Fields3d<fields_t>, typename C::ip, dim>;
    using Current_t = Current<typename C::order, dim, IP, Fields3d<fields_t>>;
    using AdvanceParticle_t = AdvanceParticle<particle_t::real_t, dim>;

    c_prm_set(prts.grid());
    real_t dqs = .5f * prts.grid().eta * prts.grid().dt;
  
    AdvanceParticle_t advance(prts.grid().dt);
    IP ip;
    Current_t c(prts.grid());

    Fields3d<fields_t> EM(flds); // FIXME, EM and J are identical here
    Fields3d<fields_t> J(flds);

    PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
      particle_t *part = &*prt_iter;
      real_t *x = &part->xi;
      real_t vv[3];

      // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt
      // FIELD INTERPOLATION

      real_t xm[3];
      for (int d = 0; d < 3; d++) {
	xm[d] = x[d] * c_prm.dxi[d];
      }
      ip.set_coeffs(xm);

      c.charge_before(ip);

      real_t E[3] = { ip.ex(EM), ip.ey(EM), ip.ez(EM) };
      real_t H[3] = { ip.hx(EM), ip.hy(EM), ip.hz(EM) };

      // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
      real_t dq = dqs * prts.prt_qni(*part) / prts.prt_mni(*part);
      advance.push_p(&part->pxi, E, H, dq);

      // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
      advance.calc_v(vv, &part->pxi);
      advance.push_x(x, vv);

      // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt
      c.charge_after(x);

      // CURRENT DENSITY AT (n+1.0)*dt
      c.prep(prts.prt_qni_wni(*part), vv);
      c.calc(J);
    }
  }

};

// ----------------------------------------------------------------------
// PscPushParticles_

template<typename PushParticles_t>
struct PscPushParticles_
{
  using Mparticles = typename PushParticles_t::Mparticles;
  using Mfields = typename PushParticles_t::Mfields;
  using mparticles_t = PscMparticles<Mparticles>;
  using mfields_t = PscMfields<Mfields>;

  static void push_mprts(struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base)
  {
    mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
    mparticles_t mp = mparticles_t(mprts);

    PushParticles_t::push_mprts(*mp.sub(), *mf.sub());

    mf.put_as(mflds_base, JXI, JXI+3);
  }
};

#ifdef CONFIG
template struct PscPushParticles_<PushParticles__<CONFIG>>;
#endif
