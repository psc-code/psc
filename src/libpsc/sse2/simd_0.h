
#include <math.h>

#include <smmintrin.h>

typedef float v4s;
typedef int v4si;

static const v4s vZero = 0.f;
static const v4s vOne = 1.f;
static const v4s vTwo = 2.f;
static const v4si iZero = 0;
static const v4s v4s_signmask = -0.f;
static inline v4s
v4s_splat(float x)
{
  return x;
}

static inline v4si
v4si_splat(int x)
{
  return x;
}

static inline v4s
v4s_load(float *p)
{
  return *p;
}

static inline void
v4s_store(float *p, v4s v)
{
  *p = v;
}

static inline v4si
v4si_load(int *p)
{
  return *p;
}

static inline void
v4si_store(int *p, v4si v)
{
  *p = v;
}

static inline void
v4s_prefetch(float *p)
{
  _mm_prefetch((const char *) p, _MM_HINT_T0);
}

static inline v4s
v4s_insert(v4s v, float s, int m)
{
  return s;
}

static inline v4si
v4si_insert(v4si v, int s, int m)
{
  return s;
}

static inline float
v4s_extract(v4s v, int m)
{
  return v;
}

static inline int
v4si_extract(v4si v, int m)
{
  return v;
}

static inline v4s
v4s_sqrt(v4s v)
{
  return sqrtf(v);
}

static inline v4s
v4s_rsqrt(v4s v)
{
  return 1.f / sqrtf(v);
}

static inline v4s
v4s_recip(v4s v)
{
#if 1
  return 1.f / v;
#else
  float x0;
  _mm_store_ss(&x0, _mm_rcp_ss(_mm_load_ss(&v)));
  return (x0 + x0) - v * (x0 * x0);
#endif
}

static inline v4si
v4s_fint(v4s v)
{
  return (int) (v + 10.f) - 10;
}

static inline v4s
v4si_to_v4s(v4si v)
{
  return (v4s) v;
}

static inline v4si
v4s_cmpge(v4s a, v4s b)
{
  return (a >= b) ? -1 : 0;
}

static inline v4si
v4s_cmpeq(v4s a, v4s b)
{
  return (a == b) ? -1 : 0;
}

static inline v4si
v4si_cmpeq(v4si a, v4si b)
{
  return (a == b) ? -1 : 0;
}

static inline v4si
v4si_cmplt(v4si a, v4si b)
{
  return (a < b) ? -1 : 0;
}

static inline v4si
v4si_cmpge(v4si a, v4si b)
{
  return (a >= b) ? -1 : 0;
}

static inline v4s
v4s_blend(v4si mask, v4s a, v4s b)
{
  return mask ? a : b;
}

static inline v4si
v4si_blend(v4si mask, v4si a, v4si b)
{
  return mask ? a : b;
}

static inline v4si
v4si_and(v4si a, v4si b)
{
  return a & b;
}

static inline v4si
v4si_or(v4si a, v4si b)
{
  return a | b;
}

static inline v4si
v4si_andnot(v4si a, v4si b)
{
  return ~a & b;
}

static inline v4si
v4si_add(v4si a, v4si b)
{
  return a + b;
}

static inline v4si
v4si_sub(v4si a, v4si b)
{
  return a - b;
}

static inline v4si
v4si_mul(v4si a, v4si b)
{
  return a * b;
}

static inline v4si
v4si_abs(v4si a)
{
  return a < 0 ? -a : a;
}

static inline v4s
v4s_floor(v4s a)
{
  return floorf(a);
}

static inline v4si
v4s_to_v4si(v4s a)
{
  return a;
}

static inline v4s
v4s_gather(real * __restrict pp, v4si off)
{
  return pp[off];
}

static inline v4si
v4si_sll(v4si a, int count)
{
  return a << count;
}

static inline v4si
v4si_cast(v4s s)
{
  union { v4s s; v4si si; } u;
  u.s = s;
  return u.si;
}

static inline v4s
v4s_cast(v4s si)
{
  union { v4s s; v4si si; } u;
  u.si = si;
  return u.s;
}

static inline v4s
v4s_xor(v4s a, v4s b)
{
  return v4s_cast(v4si_cast(a) ^ v4si_cast(b));
}

static inline v4s
v4s_and(v4s a, v4s b)
{
  return v4s_cast(v4si_cast(a) & v4si_cast(b));
}

static inline v4s
v4s_abs(v4s a)
{
  return fabsf(a);
}

