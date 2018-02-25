
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

static inline void
find_l_minmax(int *l1min, int *l1max, int k1, int lg1)
{
#if ORDER == ORDER_1ST
  if (k1 == lg1) {
    *l1min = 0; *l1max = +1;
  } else if (k1 == lg1 - 1) {
    *l1min = -1; *l1max = +1;
  } else { // (k1 == lg1 + 1)
    *l1min = 0; *l1max = +2;
  }
#elif ORDER == ORDER_2ND
  if (k1 == lg1) {
    *l1min = -1; *l1max = +1;
  } else if (k1 == lg1 - 1) {
    *l1min = -2; *l1max = +1;
  } else { // (k1 == lg1 + 1)
    *l1min = -1; *l1max = +2;
  }
#endif
}

// ----------------------------------------------------------------------
// charge density 

#if ORDER == ORDER_1ST

#define N_RHO 4
#define S_OFF 1

#elif ORDER == ORDER_2ND

#define N_RHO 5
#define S_OFF 2

#endif

#define S(s, off) s[off + S_OFF]

#define DEPOSIT(xx, k1, gx, d, dxi, s1x, lg1)		\
    int k1;						\
    gx.set(xx[d] * dxi);				\
    k1 = gx.l;						\
    set_S(s1x, k1-lg1, gx)

// ----------------------------------------------------------------------
// ZERO_S1

#define ZERO_S1 do {				\
    for (int i = -S_OFF; i < -S_OFF + N_RHO; i++) {	\
      IF_DIM_X( S(s1x, i) = 0.f; );		\
      IF_DIM_Y( S(s1y, i) = 0.f; );		\
      IF_DIM_Z( S(s1z, i) = 0.f; );		\
    }						\
  } while (0)

// ----------------------------------------------------------------------
// SUBTR_S1_S0

#define SUBTR_S1_S0 do {			\
    for (int i = -S_OFF + 1; i <= 1; i++) {	\
      IF_DIM_X( S(s1x, i) -= S(s0x, i); );	\
      IF_DIM_Y( S(s1y, i) -= S(s0y, i); );	\
      IF_DIM_Z( S(s1z, i) -= S(s0z, i); );	\
    }						\
  } while (0)

// ----------------------------------------------------------------------
// set_S

#if ORDER == ORDER_1ST
static inline void
set_S(real_t *s0, int shift, struct ip_coeff_1st<real_t> gg)
{
  S(s0, shift  ) = gg.v0;
  S(s0, shift+1) = gg.v1;
}

#elif ORDER == ORDER_2ND

static inline void
set_S(real_t *s0, int shift, struct ip_coeff_2nd<real_t> gg)
{
  // FIXME: It appears that gm/g0/g1 can be used instead of what's calculated here
  // but it needs checking.
  real_t h = gg.h;
  S(s0, shift-1) = .5f * (1.5f-std::abs(h-1.f)) * (1.5f-std::abs(h-1.f));
  S(s0, shift  ) = .75f - std::abs(h) * std::abs(h);
  S(s0, shift+1) = .5f * (1.5f-std::abs(h+1.f)) * (1.5f-std::abs(h+1.f));
}

#endif

// ======================================================================
// current

#define CURRENT_PREP_DIM(l1min, l1max, k1, cxyz, fnqx, fnqxs)	\
    int l1min, l1max; find_l_minmax(&l1min, &l1max, k1, ip.cxyz.g.l);	\
    real_t fnqx = prts.prt_qni_wni(*part) * c_prm.fnqxs;	\

#define CURRENT_PREP							\
  IF_DIM_X( CURRENT_PREP_DIM(l1min, l1max, k1, cx, fnqx, fnqxs); );	\
  IF_DIM_Y( CURRENT_PREP_DIM(l2min, l2max, k2, cy, fnqy, fnqys); );	\
  IF_DIM_Z( CURRENT_PREP_DIM(l3min, l3max, k3, cz, fnqz, fnqzs); );	\
									\
  IF_NOT_DIM_X( real_t fnqxx = vv[0] * prts.prt_qni_wni(*part) * c_prm.fnqs; ); \
  IF_NOT_DIM_Y( real_t fnqyy = vv[1] * prts.prt_qni_wni(*part) * c_prm.fnqs; ); \
  IF_NOT_DIM_Z( real_t fnqzz = vv[2] * prts.prt_qni_wni(*part) * c_prm.fnqs; )

#define CURRENT_2ND_Y						\
  real_t jyh = 0.f;					\
								\
  for (int l2 = l2min; l2 <= l2max; l2++) {			\
    real_t wx = S(s0y, l2) + .5f * S(s1y, l2);		\
    real_t wy = S(s1y, l2);				\
    real_t wz = S(s0y, l2) + .5f * S(s1y, l2);		\
    								\
    real_t jxh = fnqxx*wx;				\
    jyh -= fnqy*wy;						\
    real_t jzh = fnqzz*wz;				\
    								\
    J(JXI, 0,ip.cy.g.l+l2,0) += jxh;					\
    J(JYI, 0,ip.cy.g.l+l2,0) += jyh;				\
    J(JZI, 0,ip.cy.g.l+l2,0) += jzh;				\
  }

#define CURRENT_2ND_Z						\
  real_t jzh = 0.f;					\
  for (int l3 = l3min; l3 <= l3max; l3++) {			\
    real_t wx = S(s0z, l3) + .5f * S(s1z, l3);		\
    real_t wy = S(s0z, l3) + .5f * S(s1z, l3);		\
    real_t wz = S(s1z, l3);				\
    								\
    real_t jxh = fnqxx*wx;				\
    real_t jyh = fnqyy*wy;				\
    jzh -= fnqz*wz;						\
    								\
    J(JXI, 0,0,ip.cz.g.l+l3) += jxh;				\
    J(JYI, 0,0,ip.cz.g.l+l3) += jyh;				\
    J(JZI, 0,0,ip.cz.g.l+l3) += jzh;				\
  }

#define CURRENT_2ND_XY							\
  for (int l2 = l2min; l2 <= l2max; l2++) {				\
    real_t jxh = 0.f;						\
    for (int l1 = l1min; l1 <= l1max; l1++) {				\
      real_t wx = S(s1x, l1) * (S(s0y, l2) + .5f*S(s1y, l2));	\
      real_t wz = S(s0x, l1) * S(s0y, l2)			\
	+ .5f * S(s1x, l1) * S(s0y, l2)					\
	+ .5f * S(s0x, l1) * S(s1y, l2)					\
	+ (1.f/3.f) * S(s1x, l1) * S(s1y, l2);				\
      									\
      jxh -= fnqx*wx;							\
      J(JXI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += jxh;				\
      J(JZI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += fnqzz * wz;			\
    }									\
  }									\
  for (int l1 = l1min; l1 <= l1max; l1++) {				\
    real_t jyh = 0.f;						\
    for (int l2 = l2min; l2 <= l2max; l2++) {				\
      real_t wy = S(s1y, l2) * (S(s0x, l1) + .5f*S(s1x, l1));	\
      									\
      jyh -= fnqy*wy;							\
      J(JYI, ip.cx.g.l+l1,ip.cy.g.l+l2,0) += jyh;				\
    }									\
  }

#define CURRENT_XZ				\
  for (int l3 = l3min; l3 <= l3max; l3++) {	\
    real_t jxh = 0.f;						\
    for (int l1 = l1min; l1 < l1max; l1++) {				\
      real_t wx = S(s1x, l1) * (S(s0z, l3) + .5f*S(s1z, l3));	\
      jxh -= fnqx*wx;							\
      J(JXI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jxh;				\
    }									\
  }									\
									\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    for (int l1 = l1min; l1 <= l1max; l1++) {				\
      real_t wy = S(s0x, l1) * S(s0z, l3)			\
	+ .5f * S(s1x, l1) * S(s0z, l3)					\
	+ .5f * S(s0x, l1) * S(s1z, l3)					\
	+ (1.f/3.f) * S(s1x, l1) * S(s1z, l3);				\
      real_t jyh = fnqyy * wy;					\
      J(JYI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jyh;				\
    }									\
  }									\
  for (int l1 = l1min; l1 <= l1max; l1++) {				\
    real_t jzh = 0.f;						\
    for (int l3 = l3min; l3 < l3max; l3++) {				\
      real_t wz = S(s1z, l3) * (S(s0x, l1) + .5f*S(s1x, l1));	\
      jzh -= fnqz*wz;							\
      J(JZI, ip.cx.g.l+l1,0,ip.cz.g.l+l3) += jzh;				\
    }									\
  }

#define CURRENT_1ST_YZ							\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    for (int l2 = l2min; l2 <= l2max; l2++) {				\
      real_t wx = S(s0y, l2) * S(s0z, l3)			\
	+ .5f * S(s1y, l2) * S(s0z, l3)					\
	+ .5f * S(s0y, l2) * S(s1z, l3)					\
	+ (1.f/3.f) * S(s1y, l2) * S(s1z, l3);				\
      real_t jxh = fnqxx * wx;					\
      J(JXI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;				\
    }									\
  }									\
  									\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    real_t jyh = 0.f;						\
    for (int l2 = l2min; l2 < l2max; l2++) {				\
      real_t wy = S(s1y, l2) * (S(s0z, l3) + .5f*S(s1z, l3));	\
      jyh -= fnqy*wy;							\
      J(JYI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;				\
    }									\
  }									\
									\
  for (int l2 = l2min; l2 <= l2max; l2++) {				\
    real_t jzh = 0.f;						\
    for (int l3 = l3min; l3 < l3max; l3++) {				\
      real_t wz = S(s1z, l3) * (S(s0y, l2) + .5f*S(s1y, l2));	\
      jzh -= fnqz*wz;							\
      J(JZI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jzh;				\
    }									\
  }

#define JZH(i) jzh[i+2]
#define CURRENT_2ND_YZ							\
    real_t jxh;						\
    real_t jyh;						\
    real_t jzh[5];						\
									\
    for (int l2 = l2min; l2 <= l2max; l2++) {				\
      JZH(l2) = 0.f;							\
    }									\
    for (int l3 = l3min; l3 <= l3max; l3++) {				\
      jyh = 0.f;							\
      for (int l2 = l2min; l2 <= l2max; l2++) {				\
	real_t wx = S(s0y, l2) * S(s0z, l3)			\
	  + .5f * S(s1y, l2) * S(s0z, l3)				\
	  + .5f * S(s0y, l2) * S(s1z, l3)				\
	+ (1.f/3.f) * S(s1y, l2) * S(s1z, l3);				\
	real_t wy = S(s1y, l2) * (S(s0z, l3) + .5f*S(s1z, l3)); \
	real_t wz = S(s1z, l3) * (S(s0y, l2) + .5f*S(s1y, l2)); \
									\
	jxh = fnqxx*wx;							\
	jyh -= fnqy*wy;							\
	JZH(l2) -= fnqz*wz;						\
									\
	J(JXI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;				\
	J(JYI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;				\
	J(JZI, 0,ip.cy.g.l+l2,ip.cz.g.l+l3) += JZH(l2);			\
      }									\
    }									\

#define CURRENT_2ND_XYZ							\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    for (int l2 = l2min; l2 <= l2max; l2++) {				\
      real_t jxh = 0.f;					\
      for (int l1 = l1min; l1 <= l1max; l1++) {				\
	real_t wx = S(s1x, l1) * (S(s0y, l2) * S(s0z, l3) +	\
					   .5f * S(s1y, l2) * S(s0z, l3) + \
					   .5f * S(s0y, l2) * S(s1z, l3) + \
					   (1.f/3.f) * S(s1y, l2) * S(s1z, l3)); \
									\
	jxh -= fnqx*wx;							\
	J(JXI, ip.cx.g.l+l1,ip.cy.g.l+l2,ip.cz.g.l+l3) += jxh;			\
      }									\
    }									\
  }									\
  									\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    for (int l1 = l1min; l1 <= l1max; l1++) {				\
      real_t jyh = 0.f;					\
      for (int l2 = l2min; l2 <= l2max; l2++) {				\
	real_t wy = S(s1y, l2) * (S(s0x, l1) * S(s0z, l3) +	\
					   .5f * S(s1x, l1) * S(s0z, l3) + \
					   .5f * S(s0x, l1) * S(s1z, l3) + \
					   (1.f/3.f) * S(s1x, l1)*S(s1z, l3)); \
									\
	jyh -= fnqy*wy;							\
	J(JYI, ip.cx.g.l+l1,ip.cy.g.l+l2,ip.cz.g.l+l3) += jyh;			\
      }									\
    }									\
  }									\
									\
  for (int l2 = l2min; l2 <= l2max; l2++) {				\
    for (int l1 = l1min; l1 <= l1max; l1++) {				\
      real_t jzh = 0.f;					\
      for (int l3 = l3min; l3 <= l3max; l3++) {				\
	real_t wz = S(s1z, l3) * (S(s0x, l1) * S(s0y, l2) +	\
					   .5f * S(s1x, l1) * S(s0y, l2) +\
					   .5f * S(s0x, l1) * S(s1y, l2) +\
					   (1.f/3.f) * S(s1x, l1)*S(s1y, l2)); \
									\
	jzh -= fnqz*wz;							\
	J(JZI, ip.cx.g.l+l1,ip.cy.g.l+l2,ip.cz.g.l+l3) += jzh;			\
      }									\
    }									\
  }

#if DIM == DIM_Y
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
#if (DIM & DIM_X)
    real_t s0x[N_RHO] = {}, s1x[N_RHO];
#endif
#if (DIM & DIM_Y)
    real_t s0y[N_RHO] = {}, s1y[N_RHO];
#endif
#if (DIM & DIM_Z)
    real_t s0z[N_RHO] = {}, s1z[N_RHO];
#endif

    c_prm_set(ppsc->grid());

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

      IF_DIM_X( set_S(s0x, 0, ip.cx.g); );
      IF_DIM_Y( set_S(s0y, 0, ip.cy.g); );
      IF_DIM_Z( set_S(s0z, 0, ip.cz.g); );

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
      ZERO_S1;
      IP ip2;
      IF_DIM_X( DEPOSIT(x, k1, ip2.cx.g, 0, c_prm.dxi[0], s1x, ip.cx.g.l); );
      IF_DIM_Y( DEPOSIT(x, k2, ip2.cy.g, 1, c_prm.dxi[1], s1y, ip.cy.g.l); );
      IF_DIM_Z( DEPOSIT(x, k3, ip2.cz.g, 2, c_prm.dxi[2], s1z, ip.cz.g.l); );

      // CURRENT DENSITY AT (n+1.0)*dt

      SUBTR_S1_S0;
      CURRENT_PREP;
      CURRENT;
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

