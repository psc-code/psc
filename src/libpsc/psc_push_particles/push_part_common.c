
#define ORDER_1ST 1
#define ORDER_2ND 2

#define DIM_X 1
#define DIM_Y 2
#define DIM_Z 4
#define DIM_XY (DIM_X | DIM_Y)
#define DIM_XZ (DIM_X | DIM_Z)
#define DIM_YZ (DIM_Y | DIM_Z)
#define DIM_XYZ (DIM_X | DIM_Y | DIM_Z)

// ----------------------------------------------------------------------

#if (DIM & DIM_X)
#define IF_DIM_X(s) s do{} while(0)
#else
#define IF_DIM_X(s) do{} while(0)
#endif

#if (DIM & DIM_Y)
#define IF_DIM_Y(s) s do{} while(0)
#else
#define IF_DIM_Y(s) do{} while(0)
#endif

#if (DIM & DIM_Z)
#define IF_DIM_Z(s) s do{} while(0)
#else
#define IF_DIM_Z(s) do{} while(0)
#endif

// ----------------------------------------------------------------------
// static vars

static particle_real_t xl, yl, zl;

#if ORDER == ORDER_1ST


#elif ORDER == ORDER_2ND

#if DIM == DIM_Y
#define psc_push_particles_generic_c_push_mprts psc_push_particles_generic_c_push_mprts_y
#define PROF_NAME "genc_push_mprts_y"
#elif DIM == DIM_Z
#define psc_push_particles_generic_c_push_mprts psc_push_particles_generic_c_push_mprts_z
#define PROF_NAME "genc_push_mprts_z"
#elif DIM == DIM_XY
#define psc_push_particles_generic_c_push_mprts psc_push_particles_generic_c_push_mprts_xy
#define PROF_NAME "genc_push_mprts_xy"
#elif DIM == DIM_XZ
#define psc_push_particles_generic_c_push_mprts psc_push_particles_generic_c_push_mprts_xz
#define PROF_NAME "genc_push_mprts_xz"
#elif DIM == DIM_YZ
#define psc_push_particles_generic_c_push_mprts psc_push_particles_generic_c_push_mprts_yz
#define PROF_NAME "genc_push_mprts_yz"
#elif DIM == DIM_XYZ
#define psc_push_particles_generic_c_push_mprts psc_push_particles_generic_c_push_mprts_xyz
#define PROF_NAME "genc_push_mprts_xyz"
#endif

#endif

// ----------------------------------------------------------------------
// calc_v

static inline void
calc_v(particle_real_t *v, const particle_real_t *p)
{
  particle_real_t root = 1.f / particle_real_sqrt(1.f + sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
  for (int d = 0; d < 3; d++) {
    v[d] = p[d] * root;
  }
}

// ----------------------------------------------------------------------
// push_x

static inline void
push_x(particle_real_t *x, const particle_real_t *v)
{
  IF_DIM_X( x[0] += v[0] * xl; );
  IF_DIM_Y( x[1] += v[1] * yl; );
  IF_DIM_Z( x[2] += v[2] * zl; );
}

// ----------------------------------------------------------------------
// push_p

static inline void
push_p(particle_real_t *p, particle_real_t *E, particle_real_t *H, particle_real_t dq)
{
  particle_real_t pxm = p[0] + dq * E[0];
  particle_real_t pym = p[1] + dq * E[1];
  particle_real_t pzm = p[2] + dq * E[2];

  particle_real_t root = dq / particle_real_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
  particle_real_t taux = H[0] * root;
  particle_real_t tauy = H[1] * root;
  particle_real_t tauz = H[2] * root;
  
  particle_real_t tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
  particle_real_t pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
			 (2.f*taux*tauy+2.f*tauz)*pym + 
			 (2.f*taux*tauz-2.f*tauy)*pzm)*tau;
  particle_real_t pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
			 (1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
			 (2.f*tauy*tauz+2.f*taux)*pzm)*tau;
  particle_real_t pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
			 (2.f*tauy*tauz-2.f*taux)*pym +
			 (1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
  
  p[0] = pxp + dq * E[0];
  p[1] = pyp + dq * E[1];
  p[2] = pzp + dq * E[2];
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

// ----------------------------------------------------------------------
// interpolation

#if ORDER == ORDER_1ST

#define IP_FIELD(pf, m, gx, gy, gz)					\
  (gz##0z*(gx##0x*_F3(pf, m, l##gx##1  ,0,l##gz##3  ) +			\
	   gx##1x*_F3(pf, m, l##gx##1+1,0,l##gz##3  )) +			\
   gz##1z*(gx##0x*_F3(pf, m, l##gx##1  ,0,l##gz##3+1) +			\
	   gx##1x*_F3(pf, m, l##gx##1+1,0,l##gz##3+1)))			\

#else

#if DIM == DIM_Y
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gy##my*_F3(flds, m, 0,l##gy##2-1,0) +				\
   gy##0y*_F3(flds, m, 0,l##gy##2  ,0) +				\
   gy##1y*_F3(flds, m, 0,l##gy##2+1,0))
#elif DIM == DIM_Z
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##mz*_F3(flds, m, 0,0,l##gz##3-1) +				\
   gz##0z*_F3(flds, m, 0,0,l##gz##3  ) +				\
   gz##1z*_F3(flds, m, 0,0,l##gz##3+1))
#elif DIM == DIM_XY
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gy##my*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2-1,0) +		\
	   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2-1,0) +		\
	   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2-1,0)) +		\
   gy##0y*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2  ,0) +		\
	   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2  ,0) +		\
	   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2  ,0)) +		\
   gy##1y*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2+1,0) +		\
	   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2+1,0) +		\
	   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2+1,0)))
#elif DIM == DIM_XZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##mz*(gx##mx*_F3(flds, m, l##gx##1-1,0,l##gz##3-1) +		\
	   gx##0x*_F3(flds, m, l##gx##1  ,0,l##gz##3-1) +		\
	   gx##1x*_F3(flds, m, l##gx##1+1,0,l##gz##3-1)) +		\
   gz##0z*(gx##mx*_F3(flds, m, l##gx##1-1,0,l##gz##3  ) +		\
	   gx##0x*_F3(flds, m, l##gx##1  ,0,l##gz##3  ) +		\
	   gx##1x*_F3(flds, m, l##gx##1+1,0,l##gz##3  )) +		\
   gz##1z*(gx##mx*_F3(flds, m, l##gx##1-1,0,l##gz##3+1) +		\
	   gx##0x*_F3(flds, m, l##gx##1  ,0,l##gz##3+1) +		\
	   gx##1x*_F3(flds, m, l##gx##1+1,0,l##gz##3+1)))
#elif DIM == DIM_YZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##mz*(gy##my*_F3(flds, m, 0,l##gy##2-1,l##gz##3-1) +		\
	   gy##0y*_F3(flds, m, 0,l##gy##2  ,l##gz##3-1) +		\
	   gy##1y*_F3(flds, m, 0,l##gy##2+1,l##gz##3-1)) +		\
   gz##0z*(gy##my*_F3(flds, m, 0,l##gy##2-1,l##gz##3  ) +		\
	   gy##0y*_F3(flds, m, 0,l##gy##2  ,l##gz##3  ) +		\
	   gy##1y*_F3(flds, m, 0,l##gy##2+1,l##gz##3  )) +		\
   gz##1z*(gy##my*_F3(flds, m, 0,l##gy##2-1,l##gz##3+1) +		\
	   gy##0y*_F3(flds, m, 0,l##gy##2  ,l##gz##3+1) +		\
	   gy##1y*_F3(flds, m, 0,l##gy##2+1,l##gz##3+1)))
#elif DIM == DIM_XYZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##mz*(gy##my*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3-1) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3-1) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3-1)) + \
	   gy##0y*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3-1) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3-1) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3-1)) + \
	   gy##1y*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3-1) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3-1) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3-1))) + \
   gz##0z*(gy##my*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3  ) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3  ) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3  )) + \
	   gy##0y*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3  ) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3  ) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3  )) + \
	   gy##1y*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3  ) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3  ) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3  ))) + \
   gz##1z*(gy##my*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3+1) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3+1) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3+1)) + \
	   gy##0y*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3+1) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3+1) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3+1)) + \
	   gy##1y*(gx##mx*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3+1) + \
		   gx##0x*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3+1) + \
		   gx##1x*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3+1))))

#endif

#endif // ORDER

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
// get_fint_remainder

static inline void
get_fint_remainder(int *lg1, particle_real_t *h1, particle_real_t u)
{
  int l = particle_real_fint(u);
  *lg1 = l;
  *h1 = u-l;
}

// ----------------------------------------------------------------------
// ip_coeff_2nd

static inline void
ip_coeff_2nd(int *lg1, particle_real_t *h1,
	     particle_real_t *gmx, particle_real_t *g0x, particle_real_t *g1x,
	     particle_real_t u)
{
  int l;
  particle_real_t h;
  get_nint_remainder(&l, &h, u);
  *lg1 = l;
  *h1  = h;
  *gmx = .5f * (.5f+h)*(.5f+h);
  *g0x = .75f - h*h;
  *g1x = .5f * (.5f-h)*(.5f-h);
}

// ----------------------------------------------------------------------
// ip_coeff_1st

static inline void
ip_coeff_1st(int *lg1, particle_real_t *g0x, particle_real_t *g1x,
	     particle_real_t u)
{
  int l;
  particle_real_t h;
  get_fint_remainder(&l, &h, u);
  *lg1 = l;
  *g0x = 1.f - h;
  *g1x = h;
}

// ----------------------------------------------------------------------
// set_S_2nd
//
// FIXME: It appears that gm/g0/g1 can be used instead of what's calculated here
// but it needs checking.

static inline void
set_S_2nd(particle_real_t *s0x, int shift, particle_real_t h1)
{
  S(s0x, shift-1) = .5f*(1.5f-particle_real_abs(h1-1.f))*(1.5f-particle_real_abs(h1-1.f));
  S(s0x, shift  ) = .75f-particle_real_abs(h1)*particle_real_abs(h1);
  S(s0x, shift+1) = .5f*(1.5f-particle_real_abs(h1+1.f))*(1.5f-particle_real_abs(h1+1.f));
}

// ----------------------------------------------------------------------
// set_S_1st

static inline void
set_S_1st(particle_real_t *s0x, int shift, particle_real_t g0x, particle_real_t g1x)
{
  S(s0x, shift  ) = g0x;
  S(s0x, shift+1) = g1x;
}

// ----------------------------------------------------------------------
// find_l_minmax

static inline void
find_l_minmax(int *l1min, int *l1max, int k1, int lg1)
{
  if (k1 == lg1) {
    *l1min = -1; *l1max = +1;
  } else if (k1 == lg1 - 1) {
    *l1min = -2; *l1max = +1;
  } else { // (k1 == lg1 + 1)
    *l1min = -1; *l1max = +2;
  }
}

#if ORDER == ORDER_2ND

// ----------------------------------------------------------------------
// do_genc_push_part

static void
do_genc_push_part(int p, fields_t flds, particle_range_t prts)
{
#include "push_part_common_vars.c"
  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);

    // x^n, p^n -> x^(n+.5), p^n

    particle_real_t vv[3];
    calc_v(vv, &part->pxi);
    push_x(&part->xi, vv);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt

#if (DIM & DIM_X)
    int lg1, lh1;
    particle_real_t h1, gmx, g0x, g1x, hmx, h0x, h1x;
    ip_coeff_2nd(&lg1, &h1, &gmx, &g0x, &g1x, part->xi * dxi);
    set_S_2nd(s0x, 0, h1);

    ip_coeff_2nd(&lh1, &h1, &hmx, &h0x, &h1x, part->xi * dxi - .5f);
#endif
#if (DIM & DIM_Y)
    int lg2, lh2;
    particle_real_t h2, gmy, g0y, g1y, hmy, h0y, h1y;
    ip_coeff_2nd(&lg2, &h2, &gmy, &g0y, &g1y, part->yi * dyi);
    set_S_2nd(s0y, 0, h2);

    ip_coeff_2nd(&lh2, &h2, &hmy, &h0y, &h1y, part->yi * dyi - .5f);
#endif
#if (DIM & DIM_Z)
    int lg3, lh3;
    particle_real_t h3, gmz, g0z, g1z, hmz, h0z, h1z;
    ip_coeff_2nd(&lg3, &h3, &gmz, &g0z, &g1z, part->zi * dzi);
    set_S_2nd(s0z, 0, h3);

    ip_coeff_2nd(&lh3, &h3, &hmz, &h0z, &h1z, part->zi * dzi - .5f);
#endif

    // FIELD INTERPOLATION

    particle_real_t E[3] = { IP_FIELD(flds, EX, h, g, g),
			     IP_FIELD(flds, EY, g, h, g),
			     IP_FIELD(flds, EZ, g, g, h), };
    particle_real_t H[3] = { IP_FIELD(flds, HX, g, h, h),
			     IP_FIELD(flds, HY, h, g, h),
			     IP_FIELD(flds, HZ, h, h, g), };

     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    particle_real_t dq = dqs * part->qni / part->mni;
    push_p(&part->pxi, E, H, dq);

    calc_v(vv, &part->pxi);
    push_x(&part->xi, vv);

    for (int i = -2; i <= 2; i++) {
      IF_DIM_X( S(s1x, i) = 0.f; );
      IF_DIM_Y( S(s1y, i) = 0.f; );
      IF_DIM_Z( S(s1z, i) = 0.f; );
    }

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

#if (DIM & DIM_X)
    particle_real_t xi = part->xi + vv[0] * xl;
    int k1;
    get_nint_remainder(&k1, &h1, xi * dxi);
    set_S_2nd(s1x, k1-lg1, h1);
#endif
#if (DIM & DIM_Y)
    particle_real_t yi = part->yi + vv[1] * yl;
    int k2;
    get_nint_remainder(&k2, &h2, yi * dyi);
    set_S_2nd(s1y, k2-lg2, h2);
#endif
#if (DIM & DIM_Z)
    particle_real_t zi = part->zi + vv[2] * zl;
    int k3;
    get_nint_remainder(&k3, &h3, zi * dzi);
    set_S_2nd(s1z, k3-lg3, h3);
#endif

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      IF_DIM_X( S(s1x, i) -= S(s0x, i); );
      IF_DIM_Y( S(s1y, i) -= S(s0y, i); );
      IF_DIM_Z( S(s1z, i) -= S(s0z, i); );
    }

#if (DIM & DIM_X)
    int l1min, l1max;
    find_l_minmax(&l1min, &l1max, k1, lg1);
#endif
#if (DIM & DIM_Y)
    int l2min, l2max;
    find_l_minmax(&l2min, &l2max, k2, lg2);
#endif
#if (DIM & DIM_Z)
    int l3min, l3max;
    find_l_minmax(&l3min, &l3max, k3, lg3);
#endif

#if (DIM & DIM_X)
    particle_real_t fnqx = part->qni * part->wni * fnqxs;
#else
    particle_real_t fnqxx = vv[0] * part->qni * part->wni * fnqs;
#endif
#if (DIM & DIM_Y)
    particle_real_t fnqy = part->qni * part->wni * fnqys;
#else
    particle_real_t fnqyy = vv[1] * part->qni * part->wni * fnqs;
#endif
#if (DIM & DIM_Z)
    particle_real_t fnqz = part->qni * part->wni * fnqzs;
#else
    particle_real_t fnqzz = vv[2] * part->qni * part->wni * fnqs;
#endif

#if DIM == DIM_Y
    particle_real_t jyh = 0.f;

    for (int l2 = l2min; l2 <= l2max; l2++) {
      particle_real_t wx = S(s0y, l2) + .5f * S(s1y, l2);
      particle_real_t wy = S(s1y, l2);
      particle_real_t wz = S(s0y, l2) + .5f * S(s1y, l2);
      
      particle_real_t jxh = fnqxx*wx;
      jyh -= fnqy*wy;
      particle_real_t jzh = fnqzz*wz;

      _F3(flds, JXI, 0,lg2+l2,0) += jxh;
      _F3(flds, JYI, 0,lg2+l2,0) += jyh;
      _F3(flds, JZI, 0,lg2+l2,0) += jzh;
    }

#elif DIM == DIM_Z

    particle_real_t jzh = 0.f;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      particle_real_t wx = S(s0z, l3) + .5f * S(s1z, l3);
      particle_real_t wy = S(s0z, l3) + .5f * S(s1z, l3);
      particle_real_t wz = S(s1z, l3);
      
      particle_real_t jxh = fnqxx*wx;
      particle_real_t jyh = fnqyy*wy;
      jzh -= fnqz*wz;
      
      _F3(flds, JXI, 0,0,lg3+l3) += jxh;
      _F3(flds, JYI, 0,0,lg3+l3) += jyh;
      _F3(flds, JZI, 0,0,lg3+l3) += jzh;
    }

#elif DIM == DIM_XY

    for (int l2 = l2min; l2 <= l2max; l2++) {
      particle_real_t jxh = 0.f;
      for (int l1 = l1min; l1 <= l1max; l1++) {
	particle_real_t wx = S(s1x, l1) * (S(s0y, l2) + .5f*S(s1y, l2));
	particle_real_t wz = S(s0x, l1) * S(s0y, l2)
	  + .5f * S(s1x, l1) * S(s0y, l2)
	  + .5f * S(s0x, l1) * S(s1y, l2)
	  + (1.f/3.f) * S(s1x, l1) * S(s1y, l2);

	jxh -= fnqx*wx;
	_F3(flds, JXI, lg1+l1,lg2+l2,0) += jxh;
	_F3(flds, JZI, lg1+l1,lg2+l2,0) += fnqzz * wz;
      }
    }
    for (int l1 = l1min; l1 <= l1max; l1++) {
      particle_real_t jyh = 0.f;
      for (int l2 = l2min; l2 <= l2max; l2++) {
	particle_real_t wy = S(s1y, l2) * (S(s0x, l1) + .5f*S(s1x, l1));

	jyh -= fnqy*wy;
	_F3(flds, JYI, lg1+l1,lg2+l2,0) += jyh;
      }
    }

#elif DIM == DIM_XZ

    for (int l3 = l3min; l3 <= l3max; l3++) {
      particle_real_t jxh = 0.f;
      for (int l1 = l1min; l1 < l1max; l1++) {
	particle_real_t wx = S(s1x, l1) * (S(s0z, l3) + .5f*S(s1z, l3));
	jxh -= fnqx*wx;
	_F3(flds, JXI, lg1+l1,0,lg3+l3) += jxh;
      }
    }

    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l1 = l1min; l1 <= l1max; l1++) {
	particle_real_t wy = S(s0x, l1) * S(s0z, l3)
	  + .5f * S(s1x, l1) * S(s0z, l3)
	  + .5f * S(s0x, l1) * S(s1z, l3)
	  + (1.f/3.f) * S(s1x, l1) * S(s1z, l3);
	particle_real_t jyh = fnqyy*wy;
	_F3(flds, JYI, lg1+l1,0,lg3+l3) += jyh;
      }
    }

    for (int l1 = l1min; l1 <= l1max; l1++) {
      particle_real_t jzh = 0.f;
      for (int l3 = l3min; l3 < l3max; l3++) {
	particle_real_t wz = S(s1z, l3) * (S(s0x, l1) + .5f*S(s1x, l1));
	jzh -= fnqz*wz;
	_F3(flds, JZI, lg1+l1,0,lg3+l3) += jzh;
      }
    }

#elif DIM == DIM_YZ

    particle_real_t jxh;
    particle_real_t jyh;
    particle_real_t jzh[5];

#define JZH(i) jzh[i+2]

    for (int l2 = l2min; l2 <= l2max; l2++) {
      JZH(l2) = 0.f;
    }
    for (int l3 = l3min; l3 <= l3max; l3++) {
      jyh = 0.f;
      for (int l2 = l2min; l2 <= l2max; l2++) {
	particle_real_t wx = S(s0y, l2) * S(s0z, l3)
	  + .5f * S(s1y, l2) * S(s0z, l3)
	  + .5f * S(s0y, l2) * S(s1z, l3)
	  + (1.f/3.f) * S(s1y, l2) * S(s1z, l3);
	particle_real_t wy = S(s1y, l2) * (S(s0z, l3) + .5f*S(s1z, l3));
	particle_real_t wz = S(s1z, l3) * (S(s0y, l2) + .5f*S(s1y, l2));

	jxh = fnqxx*wx;
	jyh -= fnqy*wy;
	JZH(l2) -= fnqz*wz;

	_F3(flds, JXI, 0,lg2+l2,lg3+l3) += jxh;
	_F3(flds, JYI, 0,lg2+l2,lg3+l3) += jyh;
	_F3(flds, JZI, 0,lg2+l2,lg3+l3) += JZH(l2);
      }
    }

#elif DIM == DIM_XYZ
    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l2 = l2min; l2 <= l2max; l2++) {
	particle_real_t jxh = 0.f;
	for (int l1 = l1min; l1 <= l1max; l1++) {
	  particle_real_t wx = S(s1x, l1) * (S(s0y, l2) * S(s0z, l3) +
					     .5f * S(s1y, l2) * S(s0z, l3) +
					     .5f * S(s0y, l2) * S(s1z, l3) +
					     (1.f/3.f) * S(s1y, l2) * S(s1z, l3));

	  jxh -= fnqx*wx;
	  _F3(flds, JXI, lg1+l1,lg2+l2,lg3+l3) += jxh;
	}
      }
    }

    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l1 = l1min; l1 <= l1max; l1++) {
	particle_real_t jyh = 0.f;
	for (int l2 = l2min; l2 <= l2max; l2++) {
	  particle_real_t wy = S(s1y, l2) * (S(s0x, l1) * S(s0z, l3) +
					     .5f * S(s1x, l1) * S(s0z, l3) +
					     .5f * S(s0x, l1) * S(s1z, l3) +
					     (1.f/3.f) * S(s1x, l1)*S(s1z, l3));

	  jyh -= fnqy*wy;
	  _F3(flds, JYI, lg1+l1,lg2+l2,lg3+l3) += jyh;
	}
      }
    }

    for (int l2 = l2min; l2 <= l2max; l2++) {
      for (int l1 = l1min; l1 <= l1max; l1++) {
	particle_real_t jzh = 0.f;
	for (int l3 = l3min; l3 <= l3max; l3++) {
	  particle_real_t wz = S(s1z, l3) * (S(s0x, l1) * S(s0y, l2) +
					     .5f * S(s1x, l1) * S(s0y, l2) +
					     .5f * S(s0x, l1) * S(s1y, l2) +
					     (1.f/3.f) * S(s1x, l1)*S(s1y, l2));

	  jzh -= fnqz*wz;
	  _F3(flds, JZI, lg1+l1,lg2+l2,lg3+l3) += jzh;
	}
      }
    }

#endif
  }
}

// ----------------------------------------------------------------------
// psc_push_particles_generic_c_push_mprts

void
psc_push_particles_generic_c_push_mprts(struct psc_push_particles *push,
					struct psc_mparticles *mprts,
					struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register(PROF_NAME, 1., 0, 0);
  }
  
  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    particle_range_t prts = particle_range_mprts(mprts, p);
    fields_t_zero_range(flds, JXI, JXI + 3);
    do_genc_push_part(p, flds, prts);
  }
  prof_stop(pr);
}

#endif // ORDER_2ND

