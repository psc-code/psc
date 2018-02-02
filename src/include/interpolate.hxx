
#ifndef INTERPOLATE_HXX
#define INTERPOLATE_HXX

// ----------------------------------------------------------------------
// cuda compatibility stuff

#ifndef __CUDACC__
#define __host__
#define __device__
#endif

// ----------------------------------------------------------------------
// get_fint_remainder

template<typename R>
void __host__ __device__ get_fint_remainder(int *lg, R *h, R u)
{
#ifdef  __CUDA_ARCH__
  int l = __float2int_rd(u);
#else
  int l = fint(u);
#endif
  *lg = l;
  *h = u - l;
}

// ----------------------------------------------------------------------
// get_nint_remainder

template<typename R>
static inline void
get_nint_remainder(int *lg1, R *h1, R u)
{
  int l = nint(u);
  *lg1 = l;
  *h1 = l-u;
}

// ======================================================================
// ip_coeff

// ----------------------------------------------------------------------
// ip_coeff_1st

template<typename R>
struct ip_coeff_1st
{
  __host__ __device__ void set(R u)
  {
    R h;
    
    get_fint_remainder(&l, &h, u);
    v0 = 1.f - h;
    v1 = h;
  }

  R v0, v1;
  int l;
};

// ----------------------------------------------------------------------
// ip_coeff_2nd

template<typename R>
struct ip_coeff_2nd
{
  void set(R u)
  {
    get_nint_remainder(&l, &h, u);
    vm = .5f * (.5f+h)*(.5f+h);
    v0 = .75f - h*h;
    vp = .5f * (.5f-h)*(.5f-h);
  }

  R vm, v0, vp, h;
  int l;
};

#endif

