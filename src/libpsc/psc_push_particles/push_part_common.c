
#include "inc_defs.h"

#define CACHE_EM_J 1

#define IP_DEPOSIT

using real_t = mparticles_t::real_t;

#include "fields.hxx"
#include "inc_params.c"
#include "inc_push.c"
#include "inc_cache.c"
#include "interpolate.hxx"
using IP = InterpolateEM<Fields3d<fields_t>, opt_ip, opt_dim>;

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

#define DEPOSIT(xx, k1, gx, d, dxi, s1x, lg1)			\
  do {								\
    gx.set(xx[d] * dxi);					\
    k1 = gx.l;							\
    s1x.set(k1-lg1, gx);					\
  } while(0)

// ======================================================================
// Current

template<typename Rho1d_t, typename C>
struct Current
{
  struct Dir
  {
    Rho1d_t s0 = {}, s1;
    int lg, k;
    int lmin, lmax;
    real_t fnq;

    template<typename ip_coeff_t>
    void charge_before(ip_coeff_t g)
    {
      s0.set(0, g);
      lg = g.l;
    };
  };

  struct DirInvar
  {
    real_t fnqv;

    template<typename ip_coeff_t>
    void charge_before(ip_coeff_t g) {}
  };
  
  void charge_before(const IP& ip)
  {
    x.charge_before(ip.cx.g);
    y.charge_before(ip.cy.g);
    z.charge_before(ip.cz.g);
  }

  void charge_after(real_t xx[3])
  {
    zero_s1();

    IP ip2;
    IF_DIM_X( DEPOSIT(xx, x.k, ip2.cx.g, 0, c_prm.dxi[0], x.s1, x.lg); );
    IF_DIM_Y( DEPOSIT(xx, y.k, ip2.cy.g, 1, c_prm.dxi[1], y.s1, y.lg); );
    IF_DIM_Z( DEPOSIT(xx, z.k, ip2.cz.g, 2, c_prm.dxi[2], z.s1, z.lg); );
  }
  
  void zero_s1()
  {
    IF_DIM_X( x.s1.zero(); );
    IF_DIM_Y( y.s1.zero(); );
    IF_DIM_Z( z.s1.zero(); );
  }

  void subtr_s1_s0()
  {
    IF_DIM_X( for (int i = -x.s1.S_OFF + 1; i <= 1; i++) { x.s1[i] -= x.s0[i]; } );
    IF_DIM_Y( for (int i = -y.s1.S_OFF + 1; i <= 1; i++) { y.s1[i] -= y.s0[i]; } );
    IF_DIM_Z( for (int i = -z.s1.S_OFF + 1; i <= 1; i++) { z.s1[i] -= z.s0[i]; } );
  }

#define CURRENT_PREP_DIM(l1min, l1max, k1, cxyz, fnqx, fnqxs)	\
  find_l_minmax<typename C::order>(&l1min, &l1max, k1, ip.cxyz.g.l); \
  fnqx = qni_wni * c_prm.fnqxs;					     \

  void prep(IP& ip, real_t qni_wni, real_t vv[3])
  {
    subtr_s1_s0();
    
    IF_DIM_X( CURRENT_PREP_DIM(x.lmin, x.lmax, x.k, cx, x.fnq, fnqxs); );
    IF_DIM_Y( CURRENT_PREP_DIM(y.lmin, y.lmax, y.k, cy, y.fnq, fnqys); );
    IF_DIM_Z( CURRENT_PREP_DIM(z.lmin, z.lmax, z.k, cz, z.fnq, fnqzs); );

    IF_NOT_DIM_X( x.fnqv = vv[0] * qni_wni * c_prm.fnqs; );
    IF_NOT_DIM_Y( y.fnqv = vv[1] * qni_wni * c_prm.fnqs; );
    IF_NOT_DIM_Z( z.fnqv = vv[2] * qni_wni * c_prm.fnqs; );
  }

#if (DIM & DIM_X)
  Dir x;
#else
  DirInvar x;
#endif
#if (DIM & DIM_Y)
  Dir y;
#else
  DirInvar y;
#endif
#if (DIM & DIM_Z)
  Dir z;
#else
  DirInvar z;
#endif

#define CURRENT_2ND_Y						\
  real_t jyh = 0.f;					\
								\
  for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {			\
    real_t wx = c.y.s0[l2] + .5f * c.y.s1[l2];		\
    real_t wy = c.y.s1[l2];				\
    real_t wz = c.y.s0[l2] + .5f * c.y.s1[l2];		\
    								\
    real_t jxh = c.x.fnqv*wx;				\
    jyh -= c.y.fnq*wy;						\
    real_t jzh = c.z.fnqv*wz;				\
    								\
    J(JXI, 0,ip.cy.g.l+l2,0) += jxh;					\
    J(JYI, 0,ip.cy.g.l+l2,0) += jyh;				\
    J(JZI, 0,ip.cy.g.l+l2,0) += jzh;				\
  }

#define CURRENT_2ND_Z						\
  real_t jzh = 0.f;						\
  for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {			\
    real_t wx = c.z.s0[l3] + .5f * c.z.s1[l3];			\
    real_t wy = c.z.s0[l3] + .5f * c.z.s1[l3];			\
    real_t wz = c.z.s1[l3];					\
    								\
    real_t jxh = c.x.fnqv*wx;					\
    real_t jyh = c.y.fnqv*wy;					\
    jzh -= c.z.fnq*wz;						\
    								\
    J(JXI, 0,0,ip.cz.g.l+l3) += jxh;				\
    J(JYI, 0,0,ip.cz.g.l+l3) += jyh;				\
    J(JZI, 0,0,ip.cz.g.l+l3) += jzh;				\
  }

#define CURRENT_2ND_XY							\
  for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
    real_t jxh = 0.f;							\
    for (int l1 = c.x.lmin; l1 <= c.x.lmax; l1++) {				\
      real_t wx = c.x.s1[l1] * (c.y.s0[l2] + .5f*c.y.s1[l2]);			\
      real_t wz = c.x.s0[l1] * c.y.s0[l2]					\
	+ .5f * c.x.s1[l1] * c.y.s0[l2]					\
	+ .5f * c.x.s0[l1] * c.y.s1[l2]					\
	+ (1.f/3.f) * c.x.s1[l1] * c.y.s1[l2];				\
      									\
      jxh -= c.x.fnq*wx;							\
      J(JXI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += jxh;			\
      J(JZI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += c.z.fnqv * wz;		\
    }									\
  }									\
  for (int l1 = c.x.lmin; l1 <= c.x.lmax; l1++) {				\
    real_t jyh = 0.f;							\
    for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
      real_t wy = c.y.s1[l2] * (c.x.s0[l1] + .5f*c.x.s1[l1]);			\
      									\
      jyh -= c.y.fnq*wy;							\
      J(JYI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += jyh;			\
    }									\
  }

#define CURRENT_XZ				\
  for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {	\
    real_t jxh = 0.f;						\
    for (int l1 = c.x.lmin; l1 < c.x.lmax; l1++) {				\
      real_t wx = c.x.s1[l1] * (c.z.s0[l3] + .5f*c.z.s1[l3]);	\
      jxh -= c.x.fnq*wx;							\
      J(JXI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jxh;				\
    }									\
  }									\
									\
  for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {				\
    for (int l1 = c.x.lmin; l1 <= c.x.lmax; l1++) {				\
      real_t wy = c.x.s0[l1] * c.z.s0[l3]			\
	+ .5f * c.x.s1[l1] * c.z.s0[l3]					\
	+ .5f * c.x.s0[l1] * c.z.s1[l3]					\
	+ (1.f/3.f) * c.x.s1[l1] * c.z.s1[l3];				\
      real_t jyh = c.y.fnqv * wy;					\
      J(JYI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jyh;				\
    }									\
  }									\
  for (int l1 = c.x.lmin; l1 <= c.x.lmax; l1++) {				\
    real_t jzh = 0.f;						\
    for (int l3 = c.z.lmin; l3 < c.z.lmax; l3++) {				\
      real_t wz = c.z.s1[l3] * (c.x.s0[l1] + .5f*c.x.s1[l1]);	\
      jzh -= c.z.fnq*wz;							\
      J(JZI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jzh;				\
    }									\
  }

#define CURRENT_1ST_YZ							\
  for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {				\
    for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
      real_t wx = c.y.s0[l2] * c.z.s0[l3]			\
	+ .5f * c.y.s1[l2] * c.z.s0[l3]					\
	+ .5f * c.y.s0[l2] * c.z.s1[l3]					\
	+ (1.f/3.f) * c.y.s1[l2] * c.z.s1[l3];				\
      real_t jxh = c.x.fnqv * wx;					\
      J(JXI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;				\
    }									\
  }									\
  									\
  for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {				\
    real_t jyh = 0.f;						\
    for (int l2 = c.y.lmin; l2 < c.y.lmax; l2++) {				\
      real_t wy = c.y.s1[l2] * (c.z.s0[l3] + .5f*c.z.s1[l3]);	\
      jyh -= c.y.fnq*wy;							\
      J(JYI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;				\
    }									\
  }									\
									\
  for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
    real_t jzh = 0.f;						\
    for (int l3 = c.z.lmin; l3 < c.z.lmax; l3++) {				\
      real_t wz = c.z.s1[l3] * (c.y.s0[l2] + .5f*c.y.s1[l2]);	\
      jzh -= c.z.fnq*wz;							\
      J(JZI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jzh;				\
    }									\
  }

#define JZH(i) jzh[i+2]
#define CURRENT_2ND_YZ							\
  real_t jxh;								\
  real_t jyh;								\
  real_t jzh[5];							\
									\
    for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
      JZH(l2) = 0.f;							\
    }									\
    for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {				\
      jyh = 0.f;							\
      for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
	real_t wx = c.y.s0[l2] * c.z.s0[l3]					\
	  + .5f * c.y.s1[l2] * c.z.s0[l3]					\
	  + .5f * c.y.s0[l2] * c.z.s1[l3]					\
	+ (1.f/3.f) * c.y.s1[l2] * c.z.s1[l3];				\
	real_t wy = c.y.s1[l2] * (c.z.s0[l3] + .5f*c.z.s1[l3]);			\
	real_t wz = c.z.s1[l3] * (c.y.s0[l2] + .5f*c.y.s1[l2]);			\
									\
	jxh = c.x.fnqv*wx;							\
	jyh -= c.y.fnq*wy;							\
	JZH(l2) -= c.z.fnq*wz;						\
									\
	J(JXI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;			\
	J(JYI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;			\
	J(JZI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += JZH(l2);			\
      }									\
    }									\

#ifdef XYZ
  void calc(IP& ip, Fields3d<fields_t>& J)
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
	
	J(JXI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;
	J(JYI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;
	J(JZI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += JZH(l2);
      }									
    }									
  }
#endif
  
#define CURRENT_2ND_XYZ							\
  for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {				\
    for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
      real_t jxh = 0.f;					\
      for (int l1 = c.x.lmin; l1 <= c.x.lmax; l1++) {				\
	real_t wx = c.x.s1[l1] * (c.y.s0[l2] * c.z.s0[l3] +	\
					   .5f * c.y.s1[l2] * c.z.s0[l3] + \
					   .5f * c.y.s0[l2] * c.z.s1[l3] + \
					   (1.f/3.f) * c.y.s1[l2] * c.z.s1[l3]); \
									\
	jxh -= c.x.fnq*wx;							\
	J(JXI, ip.cx.g.l+l1,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;			\
      }									\
    }									\
  }									\
  									\
  for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {				\
    for (int l1 = c.x.lmin; l1 <= c.x.lmax; l1++) {				\
      real_t jyh = 0.f;					\
      for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
	real_t wy = c.y.s1[l2] * (c.x.s0[l1] * c.z.s0[l3] +	\
					   .5f * c.x.s1[l1] * c.z.s0[l3] + \
					   .5f * c.x.s0[l1] * c.z.s1[l3] + \
					   (1.f/3.f) * c.x.s1[l1]*c.z.s1[l3]); \
									\
	jyh -= c.y.fnq*wy;							\
	J(JYI, ip.cx.g.l+l1,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;			\
      }									\
    }									\
  }									\
									\
  for (int l2 = c.y.lmin; l2 <= c.y.lmax; l2++) {				\
    for (int l1 = c.x.lmin; l1 <= c.x.lmax; l1++) {				\
      real_t jzh = 0.f;					\
      for (int l3 = c.z.lmin; l3 <= c.z.lmax; l3++) {				\
	real_t wz = c.z.s1[l3] * (c.x.s0[l1] * c.y.s0[l2] +	\
					   .5f * c.x.s1[l1] * c.y.s0[l2] +\
					   .5f * c.x.s0[l1] * c.y.s1[l2] +\
					   (1.f/3.f) * c.x.s1[l1]*c.y.s1[l2]); \
									\
	jzh -= c.z.fnq*wz;							\
	J(JZI, ip.cx.g.l+l1,ip.cy.g.l+l2,ip.cz.g.l+l3) += jzh;			\
      }									\
    }									\
  }

};

#if DIM == DIM1
#define CURRENT do {} while(0)
#elif DIM == DIM_Y
#if ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_Y
#endif
#elif DIM == DIM_Z
#if ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_Z
#endif
#elif DIM == DIM_XY
#if ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_XY
#endif
#elif DIM == DIM_XZ
#define CURRENT CURRENT_XZ
#elif DIM == DIM_YZ

#if ORDER == ORDER_1ST
#define CURRENT CURRENT_1ST_YZ
#elif ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_YZ

#endif
#elif DIM == DIM_XYZ
#if ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_XYZ
#endif
#endif

struct CacheFields
{
#if DIM == DIM_YZ
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
#endif
};

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
  using mparticles_t = typename C::mparticles_t;
  using Mparticles = typename mparticles_t::sub_t;
  using mfields_t = typename C::mfields_t;
  using Mfields = typename mfields_t::sub_t;
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
      fields_t flds = cache.from_em(mflds[p]);
      do_push_part(flds, mprts[p]);
      cache.to_j(flds, mflds[p]);
    } 
    prof_stop(pr);
  }

private:
  static void do_push_part(fields_t flds, typename mparticles_t::patch_t& prts)
  {
    using Rho1d_t = Rho1d<typename C::order>;
    using Current_t = Current<Rho1d_t, C>;

    c_prm_set(ppsc->grid());
    Current_t c;

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
      IP ip;
      ip.set_coeffs(xm);

      c.charge_before(ip);

      real_t E[3] = { ip.ex(EM), ip.ey(EM), ip.ez(EM) };
      real_t H[3] = { ip.hx(EM), ip.hy(EM), ip.hz(EM) };

      // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
      real_t dq = c_prm.dqs * prts.prt_qni(*part) / prts.prt_mni(*part);
      push_p(&part->pxi, E, H, dq);

      // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
      calc_v(vv, &part->pxi);

      // FIXME, inelegant way of pushing full dt
      push_x(x, vv, c_prm.dt);

      // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
      c.charge_after(x);

      // CURRENT DENSITY AT (n+1.0)*dt
      c.prep(ip, prts.prt_qni_wni(*part), vv);
#ifdef XYZ
      c.calc(ip, J);
#else
      CURRENT;
#endif
    }

  }

};

// ----------------------------------------------------------------------
// PscPushParticles_

template<typename PushParticles_t>
struct PscPushParticles_
{
  using mparticles_t = typename PushParticles_t::mparticles_t;
  using mfields_t = typename PushParticles_t::mfields_t;
  using Mparticles = typename mparticles_t::sub_t;
  using Mfields = typename mfields_t::sub_t;
  
  static void push_mprts(struct psc_push_particles *push,
			 struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base)
  {
    mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
    mparticles_t mp = mparticles_t(mprts);
    
    PushParticles_t::push_mprts(*mp.sub(), *mf.sub());
    
    mf.put_as(mflds_base, JXI, JXI+3);
  }

  static void push_mprts(Mparticles& mprts, Mfields& mflds)
  {
    PushParticles_t::push_mprts(mprts, mflds);
  }
};

template struct PscPushParticles_<PushParticles__<CONFIG>>;

