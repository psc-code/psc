
#include "psc_debug.h"

// ----------------------------------------------------------------------
// interpolation

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

struct ip_coeff_1st
{
  void set(particle_real_t u)
  {
    particle_real_t h;
    
    get_fint_remainder(&l, &h, u);
    v0 = 1.f - h;
    v1 = h;
  }

  particle_real_t v0, v1;
  int l;
};

struct ip_coeff_2nd
{
  void set(particle_real_t u)
  {
    get_nint_remainder(&l, &h, u);
    vm = .5f * (.5f+h)*(.5f+h);
    v0 = .75f - h*h;
    vp = .5f * (.5f-h)*(.5f-h);
  }

  particle_real_t vm, v0, vp, h;
  int l;
};

#define DEPOSIT(xx, k1, gx, d, dxi, s1x, lg1)		\
    int k1;						\
    gx.set(xx[d] * dxi);				\
    k1 = gx.l;						\
    set_S(s1x, k1-lg1, gx)

// ----------------------------------------------------------------------
// ip_coeffs

template<typename IP_COEFF>
struct ip_coeffs_std
{
  using ip_coeff_t = IP_COEFF;
  
  void set(particle_real_t xm)
  {
    g.set(xm);
    h.set(xm - .5f);
  }
  
  ip_coeff_t g;
  ip_coeff_t h;
};

using ip_coeffs_1st = ip_coeffs_std<ip_coeff_1st>;
using ip_coeffs_2nd = ip_coeffs_std<ip_coeff_2nd>;
  
struct ip_coeffs_1st_ec
{
  using ip_coeff_t = ip_coeff_1st;
  
  void set(particle_real_t xm)
  {
    g.set(xm);
  }
  
  ip_coeff_t g;
};

// ----------------------------------------------------------------------

#if ORDER == ORDER_1ST

#else // ORDER == ORDER_2ND

#if DIM == DIM_Y
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.cy.gy.vm*EM(m, 0,ip.cy.gy.l-1,0) +				\
   ip.cy.gy.v0*EM(m, 0,ip.cy.gy.l  ,0) +				\
   ip.cy.gy.vp*EM(m, 0,ip.cy.gy.l+1,0))
#elif DIM == DIM_Z
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.cz.gz.vm*EM(m, 0,0,ip.cz.gz.l-1) +				\
   ip.cz.gz.v0*EM(m, 0,0,ip.cz.gz.l  ) +				\
   ip.cz.gz.vp*EM(m, 0,0,ip.cz.gz.l+1))
#elif DIM == DIM_XY
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.cy.gy.vm*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l-1,0) +		\
	     ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l-1,0) +		\
	     ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l-1,0)) +		\
   ip.cy.gy.v0*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l  ,0) +		\
	     ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l  ,0) +		\
	     ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l  ,0)) +		\
   ip.cy.gy.vp*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l+1,0) +		\
	     ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l+1,0) +		\
	     ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l+1,0)))
#elif DIM == DIM_XZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.cz.gz.vm*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,0,ip.cz.gz.l-1) +		\
	     ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,0,ip.cz.gz.l-1) +		\
	     ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,0,ip.cz.gz.l-1)) +		\
   ip.cz.gz.v0*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,0,ip.cz.gz.l  ) +		\
	     ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,0,ip.cz.gz.l  ) +		\
	     ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,0,ip.cz.gz.l  )) +		\
   ip.cz.gz.vp*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,0,ip.cz.gz.l+1) +		\
	     ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,0,ip.cz.gz.l+1) +		\
	     ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,0,ip.cz.gz.l+1)))
#elif DIM == DIM_YZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.cz.gz.vm*(ip.cy.gy.vm*EM(m, 0,ip.cy.gy.l-1,ip.cz.gz.l-1) +		\
	     ip.cy.gy.v0*EM(m, 0,ip.cy.gy.l  ,ip.cz.gz.l-1) +		\
	     ip.cy.gy.vp*EM(m, 0,ip.cy.gy.l+1,ip.cz.gz.l-1)) +		\
   ip.cz.gz.v0*(ip.cy.gy.vm*EM(m, 0,ip.cy.gy.l-1,ip.cz.gz.l  ) +		\
	     ip.cy.gy.v0*EM(m, 0,ip.cy.gy.l  ,ip.cz.gz.l  ) +		\
	     ip.cy.gy.vp*EM(m, 0,ip.cy.gy.l+1,ip.cz.gz.l  )) +		\
   ip.cz.gz.vp*(ip.cy.gy.vm*EM(m, 0,ip.cy.gy.l-1,ip.cz.gz.l+1) +		\
	     ip.cy.gy.v0*EM(m, 0,ip.cy.gy.l  ,ip.cz.gz.l+1) +		\
	     ip.cy.gy.vp*EM(m, 0,ip.cy.gy.l+1,ip.cz.gz.l+1)))
#elif DIM == DIM_XYZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (ip.cz.gz.vm*(ip.cy.gy.vm*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l-1,ip.cz.gz.l-1) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l-1,ip.cz.gz.l-1) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l-1,ip.cz.gz.l-1)) + \
	     ip.cy.gy.v0*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l  ,ip.cz.gz.l-1) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l  ,ip.cz.gz.l-1) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l  ,ip.cz.gz.l-1)) + \
	     ip.cy.gy.vp*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l+1,ip.cz.gz.l-1) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l+1,ip.cz.gz.l-1) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l+1,ip.cz.gz.l-1))) + \
   ip.cz.gz.v0*(ip.cy.gy.vm*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l-1,ip.cz.gz.l  ) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l-1,ip.cz.gz.l  ) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l-1,ip.cz.gz.l  )) + \
	     ip.cy.gy.v0*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l  ,ip.cz.gz.l  ) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l  ,ip.cz.gz.l  ) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l  ,ip.cz.gz.l  )) + \
	     ip.cy.gy.vp*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l+1,ip.cz.gz.l  ) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l+1,ip.cz.gz.l  ) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l+1,ip.cz.gz.l  ))) + \
   ip.cz.gz.vp*(ip.cy.gy.vm*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l-1,ip.cz.gz.l+1) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l-1,ip.cz.gz.l+1) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l-1,ip.cz.gz.l+1)) + \
	     ip.cy.gy.v0*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l  ,ip.cz.gz.l+1) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l  ,ip.cz.gz.l+1) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l  ,ip.cz.gz.l+1)) + \
	     ip.cy.gy.vp*(ip.cx.gx.vm*EM(m, ip.cx.gx.l-1,ip.cy.gy.l+1,ip.cz.gz.l+1) + \
		       ip.cx.gx.v0*EM(m, ip.cx.gx.l  ,ip.cy.gy.l+1,ip.cz.gz.l+1) + \
		       ip.cx.gx.vp*EM(m, ip.cx.gx.l+1,ip.cy.gy.l+1,ip.cz.gz.l+1))))

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

#if ORDER == ORDER_1ST
static inline void
set_S(particle_real_t *s0, int shift, struct ip_coeff_1st gg)
{
  S(s0, shift  ) = gg.v0;
  S(s0, shift+1) = gg.v1;
}

#elif ORDER == ORDER_2ND

static inline void
set_S(particle_real_t *s0, int shift, struct ip_coeff_2nd gg)
{
  // FIXME: It appears that gm/g0/g1 can be used instead of what's calculated here
  // but it needs checking.
  particle_real_t h = gg.h;
  S(s0, shift-1) = .5f*(1.5f-particle_real_abs(h-1.f))*(1.5f-particle_real_abs(h-1.f));
  S(s0, shift  ) = .75f-particle_real_abs(h)*particle_real_abs(h);
  S(s0, shift+1) = .5f*(1.5f-particle_real_abs(h+1.f))*(1.5f-particle_real_abs(h+1.f));
}

#endif

template<typename F, typename IP, typename IP_COEFFS, typename dim>
struct InterpolateEM_Helper
{
  using real_t = particle_real_t;

  static real_t ex(const IP& ip, F EM) { assert(0); }
  static real_t ey(const IP& ip, F EM) { assert(0); }
  static real_t ez(const IP& ip, F EM) { assert(0); }
  static real_t hx(const IP& ip, F EM) { assert(0); }
  static real_t hy(const IP& ip, F EM) { assert(0); }
  static real_t hz(const IP& ip, F EM) { assert(0); }
};

// ----------------------------------------------------------------------
// InterpolateEM_Helper: 1st EC, xyz

template<typename F, typename IP>
struct InterpolateEM_Helper<F, IP, ip_coeffs_1st_ec, dim_xyz>
{
  using real_t = particle_real_t;

  static real_t ex(const IP& ip, F EM)
  {
    return (ip.cz.g.v0*(ip.cy.g.v0*EM(EX, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l  ) +
			ip.cy.g.v1*EM(EX, ip.cx.g.l  ,ip.cy.g.l+1,ip.cz.g.l  )) +
	    ip.cz.g.v1*(ip.cy.g.v0*EM(EX, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l+1) +
			ip.cy.g.v1*EM(EX, ip.cx.g.l  ,ip.cy.g.l+1,ip.cz.g.l+1)));
  }

  static real_t ey(const IP& ip, F EM)
  {
    return (ip.cx.g.v0*(ip.cz.g.v0*EM(EY, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l  ) +
			ip.cz.g.v1*EM(EY, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l+1)) +
	    ip.cx.g.v1*(ip.cz.g.v0*EM(EY, ip.cx.g.l+1,ip.cy.g.l  ,ip.cz.g.l  ) +
			ip.cz.g.v1*EM(EY, ip.cx.g.l+1,ip.cy.g.l  ,ip.cz.g.l+1)));
  }

  static real_t ez(const IP& ip, F EM)
  {
    return (ip.cy.g.v0*(ip.cx.g.v0*EM(EZ, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l  ) +
			ip.cx.g.v1*EM(EZ, ip.cx.g.l+1,ip.cy.g.l  ,ip.cz.g.l  )) +
	    ip.cy.g.v1*(ip.cx.g.v0*EM(EZ, ip.cx.g.l  ,ip.cy.g.l+1,ip.cz.g.l  ) +
			ip.cx.g.v1*EM(EZ, ip.cx.g.l+1,ip.cy.g.l+1,ip.cz.g.l  )));
  }

  static real_t hx(const IP& ip, F EM)
  {
    return (ip.cx.g.v0*EM(HX, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l  ) +
	    ip.cx.g.v1*EM(HX, ip.cx.g.l+1,ip.cy.g.l  ,ip.cz.g.l  ));
  }

  static real_t hy(const IP& ip, F EM)
  {
    return (ip.cy.g.v0*EM(HY, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l  ) +
	    ip.cy.g.v1*EM(HY, ip.cx.g.l  ,ip.cy.g.l+1,ip.cz.g.l  ));	     
  }

  static real_t hz(const IP& ip, F EM)
  {
    return (ip.cz.g.v0*EM(HZ, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l  ) +
	    ip.cz.g.v1*EM(HZ, ip.cx.g.l  ,ip.cy.g.l  ,ip.cz.g.l+1));	     
  }
};

// ----------------------------------------------------------------------
// InterpolateEM_Helper: 1st EC, yz

template<typename F, typename IP>
struct InterpolateEM_Helper<F, IP, ip_coeffs_1st_ec, dim_yz>
{
  using real_t = particle_real_t;

  static real_t ex(const IP& ip, F EM)
  {
    return (ip.cz.g.v0*(ip.cy.g.v0*EM(EX, 0,ip.cy.g.l  ,ip.cz.g.l  ) +
			ip.cy.g.v1*EM(EX, 0,ip.cy.g.l+1,ip.cz.g.l  )) +
	    ip.cz.g.v1*(ip.cy.g.v0*EM(EX, 0,ip.cy.g.l  ,ip.cz.g.l+1) +
			ip.cy.g.v1*EM(EX, 0,ip.cy.g.l+1,ip.cz.g.l+1)));
  }

  static real_t ey(const IP& ip, F EM)
  {
    return (ip.cz.g.v0*EM(EY, 0,ip.cy.g.l  ,ip.cz.g.l  ) +
	    ip.cz.g.v1*EM(EY, 0,ip.cy.g.l  ,ip.cz.g.l+1));
  }

  static real_t ez(const IP& ip, F EM)
  {
    return (ip.cy.g.v0*EM(EZ, 0,ip.cy.g.l  ,ip.cz.g.l  ) +
	    ip.cy.g.v1*EM(EZ, 0,ip.cy.g.l+1,ip.cz.g.l  ));
  }

  static real_t hx(const IP& ip, F EM)
  {
    return (EM(HX, 0,ip.cy.g.l  ,ip.cz.g.l  ));
  }

  static real_t hy(const IP& ip, F EM)
  {
    return (ip.cy.g.v0*EM(HY, 0,ip.cy.g.l  ,ip.cz.g.l  ) +
	    ip.cy.g.v1*EM(HY, 0,ip.cy.g.l+1,ip.cz.g.l  ));
  }

  static real_t hz(const IP& ip, F EM)
  {
    return (ip.cz.g.v0*EM(HZ, 0,ip.cy.g.l  ,ip.cz.g.l  ) +
	    ip.cz.g.v1*EM(HZ, 0,ip.cy.g.l  ,ip.cz.g.l+1));
  }
};

// ----------------------------------------------------------------------
// InterpolateEM_Helper: 1st std, xz

template<typename F, typename IP>
struct InterpolateEM_Helper<F, IP, ip_coeffs_1st, dim_xz>
{
  using real_t = particle_real_t;
  using ip_coeff_t = typename IP::ip_coeff_t;

  static real_t cc(const ip_coeff_t& gx, const ip_coeff_t& gy, const ip_coeff_t& gz,
		   F EM, int m)
  {
    return (gz.v0*(gx.v0*EM(m, gx.l  ,0,gz.l  ) +
		   gx.v1*EM(m, gx.l+1,0,gz.l  )) +
	    gz.v1*(gx.v0*EM(m, gx.l  ,0,gz.l+1) +
		   gx.v1*EM(m, gx.l+1,0,gz.l+1)));
  }

  static real_t ex(const IP& ip, F EM) { return cc(ip.cx.h, ip.cy.g, ip.cz.g, EM, EX); }
  static real_t ey(const IP& ip, F EM) { return cc(ip.cx.g, ip.cy.h, ip.cz.g, EM, EY); }
  static real_t ez(const IP& ip, F EM) { return cc(ip.cx.g, ip.cy.g, ip.cz.h, EM, EZ); }
  static real_t hx(const IP& ip, F EM) { return cc(ip.cx.g, ip.cy.h, ip.cz.h, EM, HX); }
  static real_t hy(const IP& ip, F EM) { return cc(ip.cx.h, ip.cy.g, ip.cz.h, EM, HY); }
  static real_t hz(const IP& ip, F EM) { return cc(ip.cx.h, ip.cy.h, ip.cz.g, EM, HZ); }
};

// ----------------------------------------------------------------------
// InterpolateEM_Helper: 1st std, yz

template<typename F, typename IP>
struct InterpolateEM_Helper<F, IP, ip_coeffs_1st, dim_yz>
{
  using real_t = particle_real_t;
  using ip_coeff_t = typename IP::ip_coeff_t;

  static real_t cc(const ip_coeff_t& gx, const ip_coeff_t& gy, const ip_coeff_t& gz,
		   F EM, int m)
  {
    return (gz.v0*(gy.v0*EM(m, 0,gy.l  ,gz.l  ) +
		   gy.v1*EM(m, 0,gy.l+1,gz.l  )) +
	    gz.v1*(gy.v0*EM(m, 0,gy.l  ,gz.l+1) +
		   gy.v1*EM(m, 0,gy.l+1,gz.l+1)));
  }

  static real_t ex(const IP& ip, F EM) { return cc(ip.cx.h, ip.cy.g, ip.cz.g, EM, EX); }
  static real_t ey(const IP& ip, F EM) { return cc(ip.cx.g, ip.cy.h, ip.cz.g, EM, EY); }
  static real_t ez(const IP& ip, F EM) { return cc(ip.cx.g, ip.cy.g, ip.cz.h, EM, EZ); }
  static real_t hx(const IP& ip, F EM) { return cc(ip.cx.g, ip.cy.h, ip.cz.h, EM, HX); }
  static real_t hy(const IP& ip, F EM) { return cc(ip.cx.h, ip.cy.g, ip.cz.h, EM, HY); }
  static real_t hz(const IP& ip, F EM) { return cc(ip.cx.h, ip.cy.h, ip.cz.g, EM, HZ); }
};

// ----------------------------------------------------------------------
// dim_1, any interpolation

template<typename F, typename IP, typename IP_COEFFS>
struct InterpolateEM_Helper<F, IP, IP_COEFFS, dim_1>
{
  using real_t = particle_real_t;

  static real_t ex(const IP& ip, F EM) { return EM(EX, 0,0,0); }
  static real_t ey(const IP& ip, F EM) { return EM(EY, 0,0,0); }
  static real_t ez(const IP& ip, F EM) { return EM(EZ, 0,0,0); }
  static real_t hx(const IP& ip, F EM) { return EM(HX, 0,0,0); }
  static real_t hy(const IP& ip, F EM) { return EM(HY, 0,0,0); }
  static real_t hz(const IP& ip, F EM) { return EM(HZ, 0,0,0); }
};

// ======================================================================
// InterpolateEM

template<typename F, typename IP_COEFFS, typename dim, int N>
struct InterpolateEM
{
  using IP = InterpolateEM<F, IP_COEFFS, dim, N>;
  using ip_coeffs_t = IP_COEFFS;
  using ip_coeff_t = typename ip_coeffs_t::ip_coeff_t;
  using real_t = particle_real_t;
  
  void set_coeffs(particle_real_t xm[3])
  {
    IF_DIM_X( cx.set(xm[0]); );
    IF_DIM_Y( cy.set(xm[1]); );
    IF_DIM_Z( cz.set(xm[2]); );
  }

  using Helper = InterpolateEM_Helper<F, IP, ip_coeffs_t, dim>;
  real_t ex(F EM) { return Helper::ex(*this, EM); }
  real_t ey(F EM) { return Helper::ey(*this, EM); }
  real_t ez(F EM) { return Helper::ez(*this, EM); }
  real_t hx(F EM) { return Helper::hx(*this, EM); }
  real_t hy(F EM) { return Helper::hy(*this, EM); }
  real_t hz(F EM) { return Helper::hz(*this, EM); }
  
  particle_real_t E[3];
  particle_real_t H[3];
  ip_coeffs_t cx, cy, cz;
};

#ifndef NNN
#define NNN 0
#endif
#ifndef dim_t
#define dim_t dim_1
//#warning FIXME dim_t not defined
#endif

#if ORDER == ORDER_1ST
#if IP_VARIANT == IP_VARIANT_EC
using ip_coeffs_t = ip_coeffs_1st_ec;
#else
using ip_coeffs_t = ip_coeffs_1st;
#endif
#elif ORDER == ORDER_2ND
using ip_coeffs_t = ip_coeffs_2nd;
#endif

using IP = InterpolateEM<Fields3d<fields_t>, ip_coeffs_t, dim_t, NNN>;

#define INTERPOLATE_FIELDS(flds)					\
  ip.set_coeffs(xm);							\
  ip.E[0] = ip.ex(flds);						\
  ip.E[1] = ip.ey(flds);						\
  ip.E[2] = ip.ez(flds);						\
  ip.H[0] = ip.hx(flds);						\
  ip.H[1] = ip.hy(flds);						\
  ip.H[2] = ip.hz(flds);						\


