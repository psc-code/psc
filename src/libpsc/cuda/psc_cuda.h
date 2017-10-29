
#ifndef PSC_CUDA_H
#define PSC_CUDA_H

#include "psc.h"

#include <assert.h>
#include <math.h>
#include <psc.h>

#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

#include "cuda_bits.h"

// ======================================================================

EXTERN_C void __fields_cuda_from_device_inside(struct psc_mfields *mflds, int mb, int me);
EXTERN_C void __fields_cuda_to_device_outside(struct psc_mfields *mflds, int mb, int me);
EXTERN_C void __fields_cuda_to_device_inside(struct psc_mfields *mflds, int mb, int me);

EXTERN_C void cuda_fill_ghosts_periodic_yz(struct psc_mfields *mflds, int p, int mb, int me);
EXTERN_C void cuda_fill_ghosts_periodic_z(struct psc_mfields *mflds, int p, int mb, int me);
EXTERN_C void cuda_add_ghosts_periodic_yz(struct psc_mfields *mflds, int p, int mb, int me);
EXTERN_C void cuda_add_ghosts_periodic_z(struct psc_mfields *mflds, int p, int mb, int me);

EXTERN_C void sort_pairs_device(unsigned int *d_keys, unsigned int *d_vals, int n);
EXTERN_C void sort_pairs_host(int *d_keys, int *d_vals, int n);

#define CUDA2_STRIPE_SIZE THREADS_PER_BLOCK

EXTERN_C void psc_mparticles_cuda_get_cuda_2(struct psc_mparticles *particles,
					     void *particles_base,
					     unsigned int flags);

// These are for field boundary exchange, so they could, eventually, miss the
// interior part

#undef F3_CF_BOUNDS_CHECK

#if 0

#define F3_CF_OFF(cf, fldnr, jx,jy,jz)					\
  ((((((fldnr)								\
       * (cf)->im[2] + ((jz)-(cf)->ib[2]))				\
      * (cf)->im[1] + ((jy)-(cf)->ib[1]))				\
     * (cf)->im[0] + ((jx)-(cf)->ib[0]))))

#define F3_CF(cf, fldnr, jx,jy,jz)		\
  ((cf)->arr[F3_CF_OFF(cf, fldnr, jx,jy,jz)])

#else

#ifdef F3_CF_BOUNDS_CHECK
#define F3_CF_OFF(cf, fldnr, jx,jy,jz) ({				\
  assert(jx == 0); /* FIXME yz only! */				        \
  assert(jx >= (cf)->ib[0] && jx < (cf)->ib[0] + (cf)->im[0]);		\
  assert(jy >= (cf)->ib[1] && jy < (cf)->ib[1] + (cf)->im[1]);		\
  assert(jz >= (cf)->ib[2] && jz < (cf)->ib[2] + (cf)->im[2]);		\
  int __off = (((fldnr) * (cf)->im[2] + (jz)) * (cf)->im[1] + (jy)) * (cf)->im[0] + (jx); \
  __off; })

#else
#define F3_CF_OFF(cf, fldnr, jx,jy,jz)					\
  ((((fldnr) * (cf)->im[2] + (jz)) * (cf)->im[1] + (jy)) * (cf)->im[0] + (jx))
#endif

#define F3_CF(cf, fldnr, jx,jy,jz)					\
  ((cf)->arr_off[F3_CF_OFF(cf, fldnr, jx,jy,jz)])

#endif

#ifdef F3_CF_BOUNDS_CHECK
#define F3_CF_0_OFF(cf, fldnr, jx,jy,jz) ({				\
  assert(jx == 0); /* FIXME yz only! */				        \
  assert(jx >= 0 && jx < (cf)->im[0]);					\
  assert(jy >= 0 && jy < (cf)->im[1]);					\
  assert(jz >= 0 && jz < (cf)->im[2]);					\
  int __off = (((fldnr) * (cf)->im[2] + (jz)) * (cf)->im[1] + (jy)) * (cf)->im[0] + (jx); \
  __off; })
#else
#define F3_CF_0_OFF(cf, fldnr, jx,jy,jz)				\
  ((((fldnr) * (cf)->im[2] + (jz)) * (cf)->im[1] + (jy)) * (cf)->im[0] + (jx))
#endif

#define F3_CF_0(cf, fldnr, jx,jy,jz)					\
  ((cf)->arr[F3_CF_0_OFF(cf, fldnr, jx,jy,jz)])

#endif

