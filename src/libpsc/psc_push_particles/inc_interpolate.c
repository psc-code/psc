
#include "psc_debug.h"

// ----------------------------------------------------------------------
// interpolation

struct ip_coeff {
#if ORDER == ORDER_1ST
  particle_real_t v0, v1;
#elif ORDER == ORDER_2ND || ORDER == ORDER_1P5
  particle_real_t vm, v0, vp, h;
#endif
};

// ----------------------------------------------------------------------
// get_fint_remainder

static inline void
get_fint_remainder(int *lg, particle_real_t *h, particle_real_t u)
{
  int l = particle_real_fint(u);
  *lg = l;
  *h = u - l;
}

// ----------------------------------------------------------------------
// get_nint_remainder

static inline void
get_nint_remainder(int *lg1, particle_real_t *h1, particle_real_t u)
{
  int l = particle_real_nint(u);
  *lg1 = l;
  *h1 = l-u;
}

// ----------------------------------------------------------------------
// ip_coeff

static inline void
ip_coeff(int *lg, struct ip_coeff *gg, particle_real_t u)
{
  int l;
  particle_real_t h;

#if ORDER == ORDER_1ST
  get_fint_remainder(&l, &h, u);
  gg->v0 = 1.f - h;
  gg->v1 = h;
#elif ORDER == ORDER_2ND
  get_nint_remainder(&l, &h, u);
  gg->h  = h;
  gg->vm = .5f * (.5f+h)*(.5f+h);
  gg->v0 = .75f - h*h;
  gg->vp = .5f * (.5f-h)*(.5f-h);
#endif
  *lg = l;
}

#if IP_VARIANT != IP_VARIANT_EC
#define IP_COEFFS(lg1, lh1, gx, hx, xm) \
  ip_coeff(&lg1, &gx, xm);	\
  ip_coeff(&lh1, &hx, xm - .5f)
#else
#define IP_COEFFS(lg1, lh1, gx, hx, xm) \
  ip_coeff(&lg1, &gx, xm);
#endif
    
// ----------------------------------------------------------------------
// ip_coeff_g

static inline void
ip_coeff_g(int *lg, struct ip_coeff *gg, particle_real_t u)
{
  ip_coeff(lg, gg, u);
}

// ----------------------------------------------------------------------
// ip_coeff_h

static inline void
ip_coeff_h(int *lh, struct ip_coeff *hh, particle_real_t u)
{
#if ORDER == ORDER_1ST || ORDER == ORDER_2ND
  ip_coeff(lh, hh, u - .5f);
#elif ORDER == ORDER_1P5
  // 1+1/2 method from Sokolov paper
  // FIXME, this is almost certainly buggy
  int l;
  particle_real_t h;
  
  get_nint_remainder(&l, &h, u - .5f);

  if (h >= 0) { //0 FIXME???
    hh->vm = h;
    hh->v0 = 1.f - h;
    hh->vp = 0;
  } else { // h < 0
    hh->vm = 0;
    hh->v0 = 1.f + h;
    hh->vp = -h;
  }
  *lh = l;
  assert(0);
#endif
}

#define DEPOSIT_AND_IP_COEFFS(lg1, lh1, gx, hx, xm, s0x)	\
  ip_coeff_g(&lg1, &gx, xm);					\
  set_S(s0x, 0, gx);						\
  ip_coeff_h(&lh1, &hx, xm)
  
#define DEPOSIT(xx, k1, gx, d, dxi, s1x, lg1)		\
    int k1;						\
    ip_coeff_g(&k1, &gx, xx[d] * dxi);			\
    set_S(s1x, k1-lg1, gx)
    

// ----------------------------------------------------------------------

#if ORDER == ORDER_1ST

#if IP_VARIANT == IP_VARIANT_EC

#if DIM == DIM_1

#define IP_FIELD_EX(flds) (EM(EX, 0,0,0))
#define IP_FIELD_EY(flds) (EM(EY, 0,0,0))
#define IP_FIELD_EZ(flds) (EM(EZ, 0,0,0))
#define IP_FIELD_HX(flds) (EM(HX, 0,0,0))
#define IP_FIELD_HY(flds) (EM(HY, 0,0,0))
#define IP_FIELD_HZ(flds) (EM(HZ, 0,0,0))

#elif DIM == DIM_YZ

#define IP_FIELD_EX(flds)					\
  (ip.gz.v0*(ip.gy.v0*EM(EX, 0,ip.lg2  ,ip.lg3  ) +				\
	  ip.gy.v1*EM(EX, 0,ip.lg2+1,ip.lg3  )) +			\
   ip.gz.v1*(ip.gy.v0*EM(EX, 0,ip.lg2  ,ip.lg3+1) +				\
	  ip.gy.v1*EM(EX, 0,ip.lg2+1,ip.lg3+1)))
#define IP_FIELD_EY(flds)					\
  (ip.gz.v0*EM(EY, 0,ip.lg2  ,ip.lg3  ) +				\
   ip.gz.v1*EM(EY, 0,ip.lg2  ,ip.lg3+1))
#define IP_FIELD_EZ(flds)					\
  (ip.gy.v0*EM(EZ, 0,ip.lg2  ,ip.lg3  ) +				\
   ip.gy.v1*EM(EZ, 0,ip.lg2+1,ip.lg3  ))
#define IP_FIELD_HX(flds)			                \
  (EM(HX, 0,ip.lg2  ,ip.lg3  ))
#define IP_FIELD_HY(flds)					\
  (ip.gy.v0*EM(HY, 0,ip.lg2  ,ip.lg3  ) +				\
   ip.gy.v1*EM(HY, 0,ip.lg2+1,ip.lg3  ))
#define IP_FIELD_HZ(flds)					\
  (ip.gz.v0*EM(HZ, 0,ip.lg2  ,ip.lg3  ) +				\
   ip.gz.v1*EM(HZ, 0,ip.lg2  ,ip.lg3+1))

#elif DIM == DIM_XYZ

#define IP_FIELD_EX(flds)					\
  (ip.gz.v0*(ip.gy.v0*EM(EX, ip.lg1  ,ip.lg2  ,ip.lg3  ) +			\
	  ip.gy.v1*EM(EX, ip.lg1  ,ip.lg2+1,ip.lg3  )) +			\
   ip.gz.v1*(ip.gy.v0*EM(EX, ip.lg1  ,ip.lg2  ,ip.lg3+1) +			\
	  ip.gy.v1*EM(EX, ip.lg1  ,ip.lg2+1,ip.lg3+1)))
#define IP_FIELD_EY(flds)					\
  (ip.gx.v0*(ip.gz.v0*EM(EY, ip.lg1  ,ip.lg2  ,ip.lg3  ) +			\
	  ip.gz.v1*EM(EY, ip.lg1  ,ip.lg2  ,ip.lg3+1)) +			\
   ip.gx.v1*(ip.gz.v0*EM(EY, ip.lg1+1,ip.lg2  ,ip.lg3  ) +			\
	  ip.gz.v1*EM(EY, ip.lg1+1,ip.lg2  ,ip.lg3+1)))	     
#define IP_FIELD_EZ(flds)					\
  (ip.gy.v0*(ip.gx.v0*EM(EZ, ip.lg1  ,ip.lg2  ,ip.lg3  ) +			\
	  ip.gx.v1*EM(EZ, ip.lg1+1,ip.lg2  ,ip.lg3  )) +			\
   ip.gy.v1*(ip.gx.v0*EM(EZ, ip.lg1  ,ip.lg2+1,ip.lg3  ) +			\
	  ip.gx.v1*EM(EZ, ip.lg1+1,ip.lg2+1,ip.lg3  )))
#define IP_FIELD_HX(flds)			\
  (ip.gx.v0*EM(HX, ip.lg1  ,ip.lg2  ,ip.lg3  ) +		\
   ip.gx.v1*EM(HX, ip.lg1+1,ip.lg2  ,ip.lg3  ))
#define IP_FIELD_HY(flds)			\
  (ip.gy.v0*EM(HY, ip.lg1  ,ip.lg2  ,ip.lg3  ) +		\
   ip.gy.v1*EM(HY, ip.lg1  ,ip.lg2+1,ip.lg3  ))	     
#define IP_FIELD_HZ(flds)			\
  (ip.gz.v0*EM(HZ, ip.lg1  ,ip.lg2  ,ip.lg3  ) +		\
   ip.gz.v1*EM(HZ, ip.lg1  ,ip.lg2  ,ip.lg3+1))	     

#endif

#else // IP_VARIANT standard

#if DIM == DIM_XZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.gz##z.v0*(ip.gx##x.v0*EM(m, ip.l##gx##1  ,0,ip.l##gz##3  ) +		\
	     ip.gx##x.v1*EM(m, ip.l##gx##1+1,0,ip.l##gz##3  )) +		\
   ip.gz##z.v1*(ip.gx##x.v0*EM(m, ip.l##gx##1  ,0,ip.l##gz##3+1) +		\
	     ip.gx##x.v1*EM(m, ip.l##gx##1+1,0,ip.l##gz##3+1)))
#elif DIM == DIM_YZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.gz##z.v0*(ip.gy##y.v0*EM(m, 0,ip.l##gy##2  ,ip.l##gz##3  ) +			\
	     ip.gy##y.v1*EM(m, 0,ip.l##gy##2+1,ip.l##gz##3  )) +			\
   ip.gz##z.v1*(ip.gy##y.v0*EM(m, 0,ip.l##gy##2  ,ip.l##gz##3+1) +			\
	     ip.gy##y.v1*EM(m, 0,ip.l##gy##2+1,ip.l##gz##3+1)))
#endif

#endif // IP_VARIANT

#else // ORDER == ORDER_2ND or ORDER_1P5

#if DIM == DIM_Y
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.gy##y.vm*EM(m, 0,ip.l##gy##2-1,0) +				\
   ip.gy##y.v0*EM(m, 0,ip.l##gy##2  ,0) +				\
   ip.gy##y.vp*EM(m, 0,ip.l##gy##2+1,0))
#elif DIM == DIM_Z
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.gz##z.vm*EM(m, 0,0,ip.l##gz##3-1) +				\
   ip.gz##z.v0*EM(m, 0,0,ip.l##gz##3  ) +				\
   ip.gz##z.vp*EM(m, 0,0,ip.l##gz##3+1))
#elif DIM == DIM_XY
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.gy##y.vm*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2-1,0) +		\
	     ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2-1,0) +		\
	     ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2-1,0)) +		\
   ip.gy##y.v0*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2  ,0) +		\
	     ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2  ,0) +		\
	     ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2  ,0)) +		\
   ip.gy##y.vp*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2+1,0) +		\
	     ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2+1,0) +		\
	     ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2+1,0)))
#elif DIM == DIM_XZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.gz##z.vm*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,0,ip.l##gz##3-1) +		\
	     ip.gx##x.v0*EM(m, ip.l##gx##1  ,0,ip.l##gz##3-1) +		\
	     ip.gx##x.vp*EM(m, ip.l##gx##1+1,0,ip.l##gz##3-1)) +		\
   ip.gz##z.v0*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,0,ip.l##gz##3  ) +		\
	     ip.gx##x.v0*EM(m, ip.l##gx##1  ,0,ip.l##gz##3  ) +		\
	     ip.gx##x.vp*EM(m, ip.l##gx##1+1,0,ip.l##gz##3  )) +		\
   ip.gz##z.vp*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,0,ip.l##gz##3+1) +		\
	     ip.gx##x.v0*EM(m, ip.l##gx##1  ,0,ip.l##gz##3+1) +		\
	     ip.gx##x.vp*EM(m, ip.l##gx##1+1,0,ip.l##gz##3+1)))
#elif DIM == DIM_YZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.gz##z.vm*(ip.gy##y.vm*EM(m, 0,ip.l##gy##2-1,ip.l##gz##3-1) +		\
	     ip.gy##y.v0*EM(m, 0,ip.l##gy##2  ,ip.l##gz##3-1) +		\
	     ip.gy##y.vp*EM(m, 0,ip.l##gy##2+1,ip.l##gz##3-1)) +		\
   ip.gz##z.v0*(ip.gy##y.vm*EM(m, 0,ip.l##gy##2-1,ip.l##gz##3  ) +		\
	     ip.gy##y.v0*EM(m, 0,ip.l##gy##2  ,ip.l##gz##3  ) +		\
	     ip.gy##y.vp*EM(m, 0,ip.l##gy##2+1,ip.l##gz##3  )) +		\
   ip.gz##z.vp*(ip.gy##y.vm*EM(m, 0,ip.l##gy##2-1,ip.l##gz##3+1) +		\
	     ip.gy##y.v0*EM(m, 0,ip.l##gy##2  ,ip.l##gz##3+1) +		\
	     ip.gy##y.vp*EM(m, 0,ip.l##gy##2+1,ip.l##gz##3+1)))
#elif DIM == DIM_XYZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.gz##z.vm*(ip.gy##y.vm*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2-1,ip.l##gz##3-1) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2-1,ip.l##gz##3-1) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2-1,ip.l##gz##3-1)) + \
	     ip.gy##y.v0*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2  ,ip.l##gz##3-1) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2  ,ip.l##gz##3-1) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2  ,ip.l##gz##3-1)) + \
	     ip.gy##y.vp*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2+1,ip.l##gz##3-1) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2+1,ip.l##gz##3-1) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2+1,ip.l##gz##3-1))) + \
   ip.gz##z.v0*(ip.gy##y.vm*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2-1,ip.l##gz##3  ) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2-1,ip.l##gz##3  ) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2-1,ip.l##gz##3  )) + \
	     ip.gy##y.v0*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2  ,ip.l##gz##3  ) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2  ,ip.l##gz##3  ) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2  ,ip.l##gz##3  )) + \
	     ip.gy##y.vp*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2+1,ip.l##gz##3  ) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2+1,ip.l##gz##3  ) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2+1,ip.l##gz##3  ))) + \
   ip.gz##z.vp*(ip.gy##y.vm*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2-1,ip.l##gz##3+1) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2-1,ip.l##gz##3+1) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2-1,ip.l##gz##3+1)) + \
	     ip.gy##y.v0*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2  ,ip.l##gz##3+1) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2  ,ip.l##gz##3+1) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2  ,ip.l##gz##3+1)) + \
	     ip.gy##y.vp*(ip.gx##x.vm*EM(m, ip.l##gx##1-1,ip.l##gy##2+1,ip.l##gz##3+1) + \
		       ip.gx##x.v0*EM(m, ip.l##gx##1  ,ip.l##gy##2+1,ip.l##gz##3+1) + \
		       ip.gx##x.vp*EM(m, ip.l##gx##1+1,ip.l##gy##2+1,ip.l##gz##3+1))))

#endif

#endif // ORDER

// ======================================================================
// IP_VARIANT SFF

#if IP_VARIANT == IP_VARIANT_SFF

// FIXME, calculation of f_avg could be done at the level where we do caching, too
// (if that survives, anyway...)

#define IP_VARIANT_SFF_PREP					\
  struct psc_patch *patch = &ppsc->patch[p];			\
  								\
  /* FIXME, eventually no ghost points should be needed (?) */	\
  fields_t flds_em = fields_t((int[3]) { -1, 0, -1 },		\
			      (int[3]) { patch->ldims[0] + 2, 1, patch->ldims[2] + 1 },	\
			      6);					\
  Fields3d<fields_t> F_EM(flds_em);					\
									\
  for (int iz = -1; iz < patch->ldims[2] + 1; iz++) {			\
    for (int ix = -1; ix < patch->ldims[0] + 1; ix++) {			\
      F_EM(0, ix,0,iz) = .5 * (EM(EX, ix,0,iz) + EM(EX, ix-1,0,iz)); \
      F_EM(1, ix,0,iz) = EM(EY, ix,0,iz);		\
      F_EM(2, ix,0,iz) = .5 * (EM(EZ, ix,0,iz) + EM(EZ, ix,0,iz-1)); \
      F_EM(3, ix,0,iz) = .5 * (EM(HX, ix,0,iz) + EM(HX, ix,0,iz-1)); \
      F_EM(4, ix,0,iz) = .25 * (EM(HY, ix  ,0,iz) + EM(HY, ix  ,0,iz-1) + \
				EM(HY, ix-1,0,iz) + EM(HY, ix-1,0,iz-1)); \
      F_EM(5, ix,0,iz) = .5 * (EM(HZ, ix,0,iz) + EM(HZ, ix-1,0,iz)); \
    }									\
  }

#define IP_VARIANT_SFF_POST			\
  flds_em.dtor()

/* FIXME, we don't really need h coeffs in this case, either, though
 * the compiler may be smart enough to figure that out */

#define IP_FIELD_EX(flds) IP_FIELD(flds, EX-EX, g, g, g)
#define IP_FIELD_EY(flds) IP_FIELD(flds, EY-EX, g, g, g)
#define IP_FIELD_EZ(flds) IP_FIELD(flds, EZ-EX, g, g, g)
#define IP_FIELD_HX(flds) IP_FIELD(flds, HX-EX, g, g, g)
#define IP_FIELD_HY(flds) IP_FIELD(flds, HY-EX, g, g, g)
#define IP_FIELD_HZ(flds) IP_FIELD(flds, HZ-EX, g, g, g)

#elif IP_VARIANT == IP_VARIANT_EC
// IP_FIELD_* has already been defined earlier
#else

#define IP_VARIANT_SFF_PREP do {} while (0)
#define IP_VARIANT_SFF_POST do {} while (0)

#define IP_FIELD_EX(flds) IP_FIELD(flds, EX, h, g, g)
#define IP_FIELD_EY(flds) IP_FIELD(flds, EY, g, h, g)
#define IP_FIELD_EZ(flds) IP_FIELD(flds, EZ, g, g, h)
#define IP_FIELD_HX(flds) IP_FIELD(flds, HX, g, h, h)
#define IP_FIELD_HY(flds) IP_FIELD(flds, HY, h, g, h)
#define IP_FIELD_HZ(flds) IP_FIELD(flds, HZ, h, h, g)

#endif

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

static inline void
set_S(particle_real_t *s0, int shift, struct ip_coeff gg)
{
#if ORDER == ORDER_1ST
  S(s0, shift  ) = gg.v0;
  S(s0, shift+1) = gg.v1;
#elif ORDER == ORDER_2ND
  // FIXME: It appears that gm/g0/g1 can be used instead of what's calculated here
  // but it needs checking.
  particle_real_t h = gg.h;
  S(s0, shift-1) = .5f*(1.5f-particle_real_abs(h-1.f))*(1.5f-particle_real_abs(h-1.f));
  S(s0, shift  ) = .75f-particle_real_abs(h)*particle_real_abs(h);
  S(s0, shift+1) = .5f*(1.5f-particle_real_abs(h+1.f))*(1.5f-particle_real_abs(h+1.f));
#endif
}

#ifndef NNN
#define NNN 0
#endif

template<int N>
struct IP
{
#ifdef IP_DEPOSIT
  void set_coeffs(particle_real_t xm[3], particle_real_t s0x[],
		  particle_real_t s0y[], particle_real_t s0z[])
  {
    IF_DIM_X( DEPOSIT_AND_IP_COEFFS(lg1, lh1, gx, hx, xm[0], s0x); );
    IF_DIM_Y( DEPOSIT_AND_IP_COEFFS(lg2, lh2, gy, hy, xm[1], s0y); );
    IF_DIM_Z( DEPOSIT_AND_IP_COEFFS(lg3, lh3, gz, hz, xm[2], s0z); );
  }
#else
  void set_coeffs(particle_real_t xm[3])
  {
    IF_DIM_X( IP_COEFFS(lg1, lh1, gx, hx, xm[0]); );
    IF_DIM_Y( IP_COEFFS(lg2, lh2, gy, hy, xm[1]); );
    IF_DIM_Z( IP_COEFFS(lg3, lh3, gz, hz, xm[2]); );
  }
#endif
  
  particle_real_t E[3];
  particle_real_t H[3];
  int lg1, lh1;
  int lg2, lh2;
  int lg3, lh3;
  struct ip_coeff gx, gy, gz;
  struct ip_coeff hx, hy, hz;
};

#ifdef IP_DEPOSIT
#define INTERPOLATE_FIELDS(flds)					\
  ip.set_coeffs(xm, s0x, s0y, s0z);					\
  ip.E[0] = IP_FIELD_EX(flds);						\
  ip.E[1] = IP_FIELD_EY(flds);						\
  ip.E[2] = IP_FIELD_EZ(flds);						\
  ip.H[0] = IP_FIELD_HX(flds);						\
  ip.H[1] = IP_FIELD_HY(flds);						\
  ip.H[2] = IP_FIELD_HZ(flds);
#else
#define INTERPOLATE_FIELDS(flds)					\
  ip.set_coeffs(xm);							\
  ip.E[0] = IP_FIELD_EX(flds);						\
  ip.E[1] = IP_FIELD_EY(flds);						\
  ip.E[2] = IP_FIELD_EZ(flds);						\
  ip.H[0] = IP_FIELD_HX(flds);						\
  ip.H[1] = IP_FIELD_HY(flds);						\
  ip.H[2] = IP_FIELD_HZ(flds);
#endif

