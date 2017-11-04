
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

#define IP_COEFFS_G(lg1, gx, xm) \
  int lg1;		       \
  struct ip_coeff gx;	       \
  ip_coeff(&lg1, &gx, xm)

#if IP_VARIANT != IP_VARIANT_EC
#define IP_COEFFS_H(lh1, hx, xm) IP_COEFFS_G(lh1, hx, xm - .5f)
#else
#define IP_COEFFS_H(lh1, hx, xm) do {} while (0)
#endif

#define IP_COEFFS(lg1, lh1, gx, hx, xm) \
  IP_COEFFS_G(lg1, gx, xm);		\
  IP_COEFFS_H(lh1, hx, xm)
    
// ----------------------------------------------------------------------

#if ORDER == ORDER_1ST

#if IP_VARIANT == IP_VARIANT_EC

#if DIM == DIM_1

#define IP_FIELD_EX(flds) (F3_CACHE(flds, EX, 0,0,0))
#define IP_FIELD_EY(flds) (F3_CACHE(flds, EY, 0,0,0))
#define IP_FIELD_EZ(flds) (F3_CACHE(flds, EZ, 0,0,0))
#define IP_FIELD_HX(flds) (F3_CACHE(flds, HX, 0,0,0))
#define IP_FIELD_HY(flds) (F3_CACHE(flds, HY, 0,0,0))
#define IP_FIELD_HZ(flds) (F3_CACHE(flds, HZ, 0,0,0))

#elif DIM == DIM_YZ

#define IP_FIELD_EX(flds)						\
  (gz.v0*(gy.v0*F3_CACHE(flds, EX, 0,lg2  ,lg3  ) +			\
	  gy.v1*F3_CACHE(flds, EX, 0,lg2+1,lg3  )) +			\
   gz.v1*(gy.v0*F3_CACHE(flds, EX, 0,lg2  ,lg3+1) +			\
	  gy.v1*F3_CACHE(flds, EX, 0,lg2+1,lg3+1)))
#define IP_FIELD_EY(flds)						\
  (gz.v0*F3_CACHE(flds, EY, 0,lg2  ,lg3  ) +				\
   gz.v1*F3_CACHE(flds, EY, 0,lg2  ,lg3+1))
#define IP_FIELD_EZ(flds)						\
  (gy.v0*F3_CACHE(flds, EZ, 0,lg2  ,lg3  ) +				\
   gy.v1*F3_CACHE(flds, EZ, 0,lg2+1,lg3  ))
#define IP_FIELD_HX(flds)						\
  (F3_CACHE(flds, HX, 0,lg2  ,lg3  ))
#define IP_FIELD_HY(flds)						\
  (gy.v0*F3_CACHE(flds, HY, 0,lg2  ,lg3  ) +				\
   gy.v1*F3_CACHE(flds, HY, 0,lg2+1,lg3  ))
#define IP_FIELD_HZ(flds)						\
  (gz.v0*F3_CACHE(flds, HZ, 0,lg2  ,lg3  ) +				\
   gz.v1*F3_CACHE(flds, HZ, 0,lg2  ,lg3+1))

#elif DIM == DIM_XYZ

#define IP_FIELD_EX(flds)						\
  (gz.v0*(gy.v0*F3_CACHE(flds, EX, lg1  ,lg2  ,lg3  ) +			\
	  gy.v1*F3_CACHE(flds, EX, lg1  ,lg2+1,lg3  )) +		\
   gz.v1*(gy.v0*F3_CACHE(flds, EX, lg1  ,lg2  ,lg3+1) +			\
	  gy.v1*F3_CACHE(flds, EX, lg1  ,lg2+1,lg3+1)))
#define IP_FIELD_EY(flds)						\
  (gx.v0*(gz.v0*F3_CACHE(flds, EY, lg1  ,lg2  ,lg3  ) +			\
	  gz.v1*F3_CACHE(flds, EY, lg1  ,lg2  ,lg3+1)) +		\
   gx.v1*(gz.v0*F3_CACHE(flds, EY, lg1+1,lg2  ,lg3  ) +			\
	  gz.v1*F3_CACHE(flds, EY, lg1+1,lg2  ,lg3+1)))	     
#define IP_FIELD_EZ(flds)						\
  (gy.v0*(gx.v0*F3_CACHE(flds, EZ, lg1  ,lg2  ,lg3  ) +			\
	  gx.v1*F3_CACHE(flds, EZ, lg1+1,lg2  ,lg3  )) +		\
   gy.v1*(gx.v0*F3_CACHE(flds, EZ, lg1  ,lg2+1,lg3  ) +			\
	  gx.v1*F3_CACHE(flds, EZ, lg1+1,lg2+1,lg3  )))
#define IP_FIELD_HX(flds)					\
  (gx.v0*F3_CACHE(flds, HX, lg1  ,lg2  ,lg3  ) +		\
   gx.v1*F3_CACHE(flds, HX, lg1+1,lg2  ,lg3  ))
#define IP_FIELD_HY(flds)					\
  (gy.v0*F3_CACHE(flds, HY, lg1  ,lg2  ,lg3  ) +		\
   gy.v1*F3_CACHE(flds, HY, lg1  ,lg2+1,lg3  ))	     
#define IP_FIELD_HZ(flds)					\
  (gz.v0*F3_CACHE(flds, HZ, lg1  ,lg2  ,lg3  ) +		\
   gz.v1*F3_CACHE(flds, HZ, lg1  ,lg2  ,lg3+1))	     

#endif

#else // IP_VARIANT standard

#if DIM == DIM_XZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.v0*(gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3  ) +		\
	     gx##x.v1*_F3(flds, m, l##gx##1+1,0,l##gz##3  )) +		\
   gz##z.v1*(gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3+1) +		\
	     gx##x.v1*_F3(flds, m, l##gx##1+1,0,l##gz##3+1)))
#elif DIM == DIM_YZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.v0*(gy##y.v0*F3_CACHE(flds, m, 0,l##gy##2  ,l##gz##3  ) +	\
	     gy##y.v1*F3_CACHE(flds, m, 0,l##gy##2+1,l##gz##3  )) +	\
   gz##z.v1*(gy##y.v0*F3_CACHE(flds, m, 0,l##gy##2  ,l##gz##3+1) +	\
	     gy##y.v1*F3_CACHE(flds, m, 0,l##gy##2+1,l##gz##3+1)))
#endif

#endif // IP_VARIANT

#else // ORDER == ORDER_2ND or ORDER_1P5

#if DIM == DIM_Y
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gy##y.vm*_F3(flds, m, 0,l##gy##2-1,0) +				\
   gy##y.v0*_F3(flds, m, 0,l##gy##2  ,0) +				\
   gy##y.vp*_F3(flds, m, 0,l##gy##2+1,0))
#elif DIM == DIM_Z
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.vm*_F3(flds, m, 0,0,l##gz##3-1) +				\
   gz##z.v0*_F3(flds, m, 0,0,l##gz##3  ) +				\
   gz##z.vp*_F3(flds, m, 0,0,l##gz##3+1))
#elif DIM == DIM_XY
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gy##y.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2-1,0) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2-1,0) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2-1,0)) +		\
   gy##y.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2  ,0) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2  ,0) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2  ,0)) +		\
   gy##y.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2+1,0) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2+1,0) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2+1,0)))
#elif DIM == DIM_XZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,0,l##gz##3-1) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3-1) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,0,l##gz##3-1)) +		\
   gz##z.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,0,l##gz##3  ) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3  ) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,0,l##gz##3  )) +		\
   gz##z.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,0,l##gz##3+1) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3+1) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,0,l##gz##3+1)))
#elif DIM == DIM_YZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.vm*(gy##y.vm*_F3(flds, m, 0,l##gy##2-1,l##gz##3-1) +		\
	     gy##y.v0*_F3(flds, m, 0,l##gy##2  ,l##gz##3-1) +		\
	     gy##y.vp*_F3(flds, m, 0,l##gy##2+1,l##gz##3-1)) +		\
   gz##z.v0*(gy##y.vm*_F3(flds, m, 0,l##gy##2-1,l##gz##3  ) +		\
	     gy##y.v0*_F3(flds, m, 0,l##gy##2  ,l##gz##3  ) +		\
	     gy##y.vp*_F3(flds, m, 0,l##gy##2+1,l##gz##3  )) +		\
   gz##z.vp*(gy##y.vm*_F3(flds, m, 0,l##gy##2-1,l##gz##3+1) +		\
	     gy##y.v0*_F3(flds, m, 0,l##gy##2  ,l##gz##3+1) +		\
	     gy##y.vp*_F3(flds, m, 0,l##gy##2+1,l##gz##3+1)))
#elif DIM == DIM_XYZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.vm*(gy##y.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3-1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3-1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3-1)) + \
	     gy##y.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3-1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3-1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3-1)) + \
	     gy##y.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3-1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3-1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3-1))) + \
   gz##z.v0*(gy##y.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3  ) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3  ) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3  )) + \
	     gy##y.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3  ) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3  ) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3  )) + \
	     gy##y.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3  ) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3  ) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3  ))) + \
   gz##z.vp*(gy##y.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3+1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3+1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3+1)) + \
	     gy##y.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3+1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3+1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3+1)) + \
	     gy##y.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3+1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3+1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3+1))))

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
  fields_t flds_em = fields_t_ctor((int[3]) { -1, 0, -1 },		\
				    (int[3]) { patch->ldims[0] + 2, 1, patch->ldims[2] + 1 },\
				    6);					\
									\
  for (int iz = -1; iz < patch->ldims[2] + 1; iz++) {			\
    for (int ix = -1; ix < patch->ldims[0] + 1; ix++) {			\
      _F3(flds_em, 0, ix,0,iz) = .5 * (_F3(flds, EX, ix,0,iz) + _F3(flds, EX, ix-1,0,iz)); \
      _F3(flds_em, 1, ix,0,iz) = _F3(flds, EY, ix,0,iz);		\
      _F3(flds_em, 2, ix,0,iz) = .5 * (_F3(flds, EZ, ix,0,iz) + _F3(flds, EZ, ix,0,iz-1)); \
      _F3(flds_em, 3, ix,0,iz) = .5 * (_F3(flds, HX, ix,0,iz) + _F3(flds, HX, ix,0,iz-1)); \
      _F3(flds_em, 4, ix,0,iz) = .25 * (_F3(flds, HY, ix  ,0,iz) + _F3(flds, HY, ix  ,0,iz-1) + \
					 _F3(flds, HY, ix-1,0,iz) + _F3(flds, HY, ix-1,0,iz-1)); \
      _F3(flds_em, 5, ix,0,iz) = .5 * (_F3(flds, HZ, ix,0,iz) + _F3(flds, HZ, ix-1,0,iz)); \
    }									\
  }

#define IP_VARIANT_SFF_POST			\
  fields_t_dtor(&flds_em)

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

#define IP_VARIANT_SFF_PREP \
  fields_t flds_em = flds
#define IP_VARIANT_SFF_POST do {} while (0)

#define IP_FIELD_EX(flds) IP_FIELD(flds, EX, h, g, g)
#define IP_FIELD_EY(flds) IP_FIELD(flds, EY, g, h, g)
#define IP_FIELD_EZ(flds) IP_FIELD(flds, EZ, g, g, h)
#define IP_FIELD_HX(flds) IP_FIELD(flds, HX, g, h, h)
#define IP_FIELD_HY(flds) IP_FIELD(flds, HY, h, g, h)
#define IP_FIELD_HZ(flds) IP_FIELD(flds, HZ, h, h, g)

#endif

#define INTERPOLATE_FIELDS(flds)					\
  particle_real_t E[3] = { IP_FIELD_EX(flds),				\
                           IP_FIELD_EY(flds),				\
                           IP_FIELD_EZ(flds), };			\
  particle_real_t H[3] = { IP_FIELD_HX(flds),				\
                           IP_FIELD_HY(flds),				\
                           IP_FIELD_HZ(flds), }

