
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

#define DEPOSIT(xx, k1, gx, d, dxi, s1x, lg1)		\
    gx.set(xx[d] * dxi);				\
    c.k1 = gx.l;						\
    s1x.set(c.k1-lg1, gx)

// ======================================================================
// Current

template<typename Rho1d_t>
struct Current
{
  void zero_s1()
  {
    IF_DIM_X( s1x.zero(); );
    IF_DIM_Y( s1y.zero(); );
    IF_DIM_Z( s1z.zero(); );
  }

  void subtr_s1_s0()
  {
    IF_DIM_X( for (int i = -s1x.S_OFF + 1; i <= 1; i++) { s1x[i] -= s0x[i]; } );
    IF_DIM_Y( for (int i = -s1y.S_OFF + 1; i <= 1; i++) { s1y[i] -= s0y[i]; } );
    IF_DIM_Z( for (int i = -s1z.S_OFF + 1; i <= 1; i++) { s1z[i] -= s0z[i]; } );
  }

#if (DIM & DIM_X)
  Rho1d_t s0x = {}, s1x;
  int k1, l1min, l1max;
  real_t fnqx;
#else
  real_t fnqxx;
#endif
#if (DIM & DIM_Y)
  Rho1d_t s0y = {}, s1y;
  int k2, l2min, l2max;
  real_t fnqy;
#else
  real_t fnqyy;
#endif
#if (DIM & DIM_Z)
  Rho1d_t s0z = {}, s1z;
  int k3, l3min, l3max;
  real_t fnqz;
#else
  real_t fnqzz;
#endif

#define CURRENT_PREP_DIM(l1min, l1max, k1, cxyz, fnqx, fnqxs)	\
  find_l_minmax<typename C::order>(&c.l1min, &c.l1max, c.k1, ip.cxyz.g.l); \
  c.fnqx = prts.prt_qni_wni(*part) * c_prm.fnqxs;				\

#define CURRENT_PREP							\
  IF_DIM_X( CURRENT_PREP_DIM(l1min, l1max, k1, cx, fnqx, fnqxs); );	\
  IF_DIM_Y( CURRENT_PREP_DIM(l2min, l2max, k2, cy, fnqy, fnqys); );	\
  IF_DIM_Z( CURRENT_PREP_DIM(l3min, l3max, k3, cz, fnqz, fnqzs); );	\
									\
  IF_NOT_DIM_X( c.fnqxx = vv[0] * prts.prt_qni_wni(*part) * c_prm.fnqs; ); \
  IF_NOT_DIM_Y( c.fnqyy = vv[1] * prts.prt_qni_wni(*part) * c_prm.fnqs; ); \
  IF_NOT_DIM_Z( c.fnqzz = vv[2] * prts.prt_qni_wni(*part) * c_prm.fnqs; )

#define CURRENT_2ND_Y						\
  real_t jyh = 0.f;					\
								\
  for (int l2 = c.l2min; l2 <= c.l2max; l2++) {			\
    real_t wx = c.s0y[l2] + .5f * c.s1y[l2];		\
    real_t wy = c.s1y[l2];				\
    real_t wz = c.s0y[l2] + .5f * c.s1y[l2];		\
    								\
    real_t jxh = c.fnqxx*wx;				\
    jyh -= c.fnqy*wy;						\
    real_t jzh = c.fnqzz*wz;				\
    								\
    J(JXI, 0,ip.cy.g.l+l2,0) += jxh;					\
    J(JYI, 0,ip.cy.g.l+l2,0) += jyh;				\
    J(JZI, 0,ip.cy.g.l+l2,0) += jzh;				\
  }

#define CURRENT_2ND_Z						\
  real_t jzh = 0.f;						\
  for (int l3 = c.l3min; l3 <= c.l3max; l3++) {			\
    real_t wx = c.s0z[l3] + .5f * c.s1z[l3];			\
    real_t wy = c.s0z[l3] + .5f * c.s1z[l3];			\
    real_t wz = c.s1z[l3];					\
    								\
    real_t jxh = c.fnqxx*wx;					\
    real_t jyh = c.fnqyy*wy;					\
    jzh -= c.fnqz*wz;						\
    								\
    J(JXI, 0,0,ip.cz.g.l+l3) += jxh;				\
    J(JYI, 0,0,ip.cz.g.l+l3) += jyh;				\
    J(JZI, 0,0,ip.cz.g.l+l3) += jzh;				\
  }

#define CURRENT_2ND_XY							\
  for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
    real_t jxh = 0.f;							\
    for (int l1 = c.l1min; l1 <= c.l1max; l1++) {				\
      real_t wx = c.s1x[l1] * (c.s0y[l2] + .5f*c.s1y[l2]);			\
      real_t wz = c.s0x[l1] * c.s0y[l2]					\
	+ .5f * c.s1x[l1] * c.s0y[l2]					\
	+ .5f * c.s0x[l1] * c.s1y[l2]					\
	+ (1.f/3.f) * c.s1x[l1] * c.s1y[l2];				\
      									\
      jxh -= c.fnqx*wx;							\
      J(JXI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += jxh;			\
      J(JZI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += c.fnqzz * wz;		\
    }									\
  }									\
  for (int l1 = c.l1min; l1 <= c.l1max; l1++) {				\
    real_t jyh = 0.f;							\
    for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
      real_t wy = c.s1y[l2] * (c.s0x[l1] + .5f*c.s1x[l1]);			\
      									\
      jyh -= c.fnqy*wy;							\
      J(JYI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += jyh;			\
    }									\
  }

#define CURRENT_XZ				\
  for (int l3 = c.l3min; l3 <= c.l3max; l3++) {	\
    real_t jxh = 0.f;						\
    for (int l1 = c.l1min; l1 < c.l1max; l1++) {				\
      real_t wx = c.s1x[l1] * (c.s0z[l3] + .5f*c.s1z[l3]);	\
      jxh -= c.fnqx*wx;							\
      J(JXI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jxh;				\
    }									\
  }									\
									\
  for (int l3 = c.l3min; l3 <= c.l3max; l3++) {				\
    for (int l1 = c.l1min; l1 <= c.l1max; l1++) {				\
      real_t wy = c.s0x[l1] * c.s0z[l3]			\
	+ .5f * c.s1x[l1] * c.s0z[l3]					\
	+ .5f * c.s0x[l1] * c.s1z[l3]					\
	+ (1.f/3.f) * c.s1x[l1] * c.s1z[l3];				\
      real_t jyh = c.fnqyy * wy;					\
      J(JYI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jyh;				\
    }									\
  }									\
  for (int l1 = c.l1min; l1 <= c.l1max; l1++) {				\
    real_t jzh = 0.f;						\
    for (int l3 = c.l3min; l3 < c.l3max; l3++) {				\
      real_t wz = c.s1z[l3] * (c.s0x[l1] + .5f*c.s1x[l1]);	\
      jzh -= c.fnqz*wz;							\
      J(JZI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jzh;				\
    }									\
  }

#define CURRENT_1ST_YZ							\
  for (int l3 = c.l3min; l3 <= c.l3max; l3++) {				\
    for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
      real_t wx = c.s0y[l2] * c.s0z[l3]			\
	+ .5f * c.s1y[l2] * c.s0z[l3]					\
	+ .5f * c.s0y[l2] * c.s1z[l3]					\
	+ (1.f/3.f) * c.s1y[l2] * c.s1z[l3];				\
      real_t jxh = c.fnqxx * wx;					\
      J(JXI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;				\
    }									\
  }									\
  									\
  for (int l3 = c.l3min; l3 <= c.l3max; l3++) {				\
    real_t jyh = 0.f;						\
    for (int l2 = c.l2min; l2 < c.l2max; l2++) {				\
      real_t wy = c.s1y[l2] * (c.s0z[l3] + .5f*c.s1z[l3]);	\
      jyh -= c.fnqy*wy;							\
      J(JYI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;				\
    }									\
  }									\
									\
  for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
    real_t jzh = 0.f;						\
    for (int l3 = c.l3min; l3 < c.l3max; l3++) {				\
      real_t wz = c.s1z[l3] * (c.s0y[l2] + .5f*c.s1y[l2]);	\
      jzh -= c.fnqz*wz;							\
      J(JZI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jzh;				\
    }									\
  }

#define JZH(i) jzh[i+2]
#define CURRENT_2ND_YZ							\
  real_t jxh;								\
  real_t jyh;								\
  real_t jzh[5];							\
									\
    for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
      JZH(l2) = 0.f;							\
    }									\
    for (int l3 = c.l3min; l3 <= c.l3max; l3++) {				\
      jyh = 0.f;							\
      for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
	real_t wx = c.s0y[l2] * c.s0z[l3]					\
	  + .5f * c.s1y[l2] * c.s0z[l3]					\
	  + .5f * c.s0y[l2] * c.s1z[l3]					\
	+ (1.f/3.f) * c.s1y[l2] * c.s1z[l3];				\
	real_t wy = c.s1y[l2] * (c.s0z[l3] + .5f*c.s1z[l3]);			\
	real_t wz = c.s1z[l3] * (c.s0y[l2] + .5f*c.s1y[l2]);			\
									\
	jxh = c.fnqxx*wx;							\
	jyh -= c.fnqy*wy;							\
	JZH(l2) -= c.fnqz*wz;						\
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
    
    for (int l2 = l2min; l2 <= l2max; l2++) {
      JZH(l2) = 0.f;
    }
    for (int l3 = l3min; l3 <= l3max; l3++) {
      jyh = 0.f;
      for (int l2 = l2min; l2 <= l2max; l2++) {
	real_t wx = s0y[l2] * s0z[l3]
	  + .5f * s1y[l2] * s0z[l3]
	  + .5f * s0y[l2] * s1z[l3]
	  + (1.f/3.f) * s1y[l2] * s1z[l3];
	real_t wy = s1y[l2] * (s0z[l3] + .5f*s1z[l3]);
	real_t wz = s1z[l3] * (s0y[l2] + .5f*s1y[l2]);
	
	jxh = fnqxx*wx;
	jyh -= fnqy*wy;
	JZH(l2) -= fnqz*wz;
	
	J(JXI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;
	J(JYI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;
	J(JZI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += JZH(l2);
      }									
    }									
  }
#endif
  
#define CURRENT_2ND_XYZ							\
  for (int l3 = c.l3min; l3 <= c.l3max; l3++) {				\
    for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
      real_t jxh = 0.f;					\
      for (int l1 = c.l1min; l1 <= c.l1max; l1++) {				\
	real_t wx = c.s1x[l1] * (c.s0y[l2] * c.s0z[l3] +	\
					   .5f * c.s1y[l2] * c.s0z[l3] + \
					   .5f * c.s0y[l2] * c.s1z[l3] + \
					   (1.f/3.f) * c.s1y[l2] * c.s1z[l3]); \
									\
	jxh -= c.fnqx*wx;							\
	J(JXI, ip.cx.g.l+l1,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;			\
      }									\
    }									\
  }									\
  									\
  for (int l3 = c.l3min; l3 <= c.l3max; l3++) {				\
    for (int l1 = c.l1min; l1 <= c.l1max; l1++) {				\
      real_t jyh = 0.f;					\
      for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
	real_t wy = c.s1y[l2] * (c.s0x[l1] * c.s0z[l3] +	\
					   .5f * c.s1x[l1] * c.s0z[l3] + \
					   .5f * c.s0x[l1] * c.s1z[l3] + \
					   (1.f/3.f) * c.s1x[l1]*c.s1z[l3]); \
									\
	jyh -= c.fnqy*wy;							\
	J(JYI, ip.cx.g.l+l1,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;			\
      }									\
    }									\
  }									\
									\
  for (int l2 = c.l2min; l2 <= c.l2max; l2++) {				\
    for (int l1 = c.l1min; l1 <= c.l1max; l1++) {				\
      real_t jzh = 0.f;					\
      for (int l3 = c.l3min; l3 <= c.l3max; l3++) {				\
	real_t wz = c.s1z[l3] * (c.s0x[l1] * c.s0y[l2] +	\
					   .5f * c.s1x[l1] * c.s0y[l2] +\
					   .5f * c.s0x[l1] * c.s1y[l2] +\
					   (1.f/3.f) * c.s1x[l1]*c.s1y[l2]); \
									\
	jzh -= c.fnqz*wz;							\
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
    using Current_t = Current<Rho1d_t>;

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

      IF_DIM_X( c.s0x.set(0, ip.cx.g); );
      IF_DIM_Y( c.s0y.set(0, ip.cy.g); );
      IF_DIM_Z( c.s0z.set(0, ip.cz.g); );

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
      c.zero_s1();
      IP ip2;
      IF_DIM_X( DEPOSIT(x, k1, ip2.cx.g, 0, c_prm.dxi[0], c.s1x, ip.cx.g.l); );
      IF_DIM_Y( DEPOSIT(x, k2, ip2.cy.g, 1, c_prm.dxi[1], c.s1y, ip.cy.g.l); );
      IF_DIM_Z( DEPOSIT(x, k3, ip2.cz.g, 2, c_prm.dxi[2], c.s1z, ip.cz.g.l); );

      // CURRENT DENSITY AT (n+1.0)*dt

      c.subtr_s1_s0();
      CURRENT_PREP;
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

