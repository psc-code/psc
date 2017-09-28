
#ifndef PSC_CUDA_H
#define PSC_CUDA_H

#include "psc.h"

#include <assert.h>
#include <math.h>
#include <psc.h>

#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

// ======================================================================

EXTERN_C void cuda_axpy_comp_yz(struct psc_mfields *y, int ym, float a, struct psc_mfields *x, int xm, int p);
EXTERN_C void cuda_zero_comp_yz(struct psc_mfields *x, int xm, int p);

EXTERN_C void cuda_marder_correct_yz(struct psc_mfields *mflds, struct psc_mfields *mf,
				     int p, int ldims[3], float fac[3],
				     int ly[3], int ry[3],
				     int lz[3], int rz[3]);
EXTERN_C void cuda_calc_dive_yz(struct psc_mfields *mflds, struct psc_mfields *mf, int p);

EXTERN_C void yz_moments_rho_1st_nc_cuda_run_patches(struct psc_mparticles *mprts, struct psc_mfields *mres);
EXTERN_C void yz_moments_n_1st_cuda_run_patches(struct psc_mparticles *mprts, struct psc_mfields *mres);

struct d_particle {
  real xi[3];
  real kind_as_float;
  real pxi[3];
  real qni_wni;
};

#define THREADS_PER_BLOCK 256

#define LOAD_PARTICLE(pp, d_p, n) do {					\
    (pp).xi[0]         = d_p.xi4[n].x;					\
    (pp).xi[1]         = d_p.xi4[n].y;					\
    (pp).xi[2]         = d_p.xi4[n].z;					\
    (pp).kind_as_float = d_p.xi4[n].w;					\
    (pp).pxi[0]        = d_p.pxi4[n].x;					\
    (pp).pxi[1]        = d_p.pxi4[n].y;					\
    (pp).pxi[2]        = d_p.pxi4[n].z;					\
    (pp).qni_wni       = d_p.pxi4[n].w;					\
} while (0)

#define STORE_PARTICLE_POS(pp, d_p, n) do {				\
    d_p.xi4[n].x = (pp).xi[0];						\
    d_p.xi4[n].y = (pp).xi[1];						\
    d_p.xi4[n].z = (pp).xi[2];						\
    d_p.xi4[n].w = (pp).kind_as_float;					\
} while (0)

#define STORE_PARTICLE_MOM(pp, d_p, n) do {				\
    d_p.pxi4[n].x = (pp).pxi[0];					\
    d_p.pxi4[n].y = (pp).pxi[1];					\
    d_p.pxi4[n].z = (pp).pxi[2];					\
    d_p.pxi4[n].w = (pp).qni_wni;					\
} while (0)

EXTERN_C void __psc_mfields_cuda_setup(struct psc_mfields *mflds);
EXTERN_C void __psc_mfields_cuda_destroy(struct psc_mfields *mflds);
EXTERN_C void __fields_cuda_to_device(struct psc_mfields *mflds, int p, real *h_flds, int mb, int me);
EXTERN_C void __fields_cuda_from_device(struct psc_mfields *mflds, int p, real *h_flds, int mb, int me);

EXTERN_C void __fields_cuda_from_device_inside(struct psc_mfields *mflds, int mb, int me);
EXTERN_C void __fields_cuda_to_device_outside(struct psc_mfields *mflds, int mb, int me);
EXTERN_C void __fields_cuda_to_device_inside(struct psc_mfields *mflds, int mb, int me);

EXTERN_C void cuda_fill_ghosts_periodic_yz(struct psc_mfields *mflds, int p, int mb, int me);
EXTERN_C void cuda_fill_ghosts_periodic_z(struct psc_mfields *mflds, int p, int mb, int me);
EXTERN_C void cuda_add_ghosts_periodic_yz(struct psc_mfields *mflds, int p, int mb, int me);
EXTERN_C void cuda_add_ghosts_periodic_z(struct psc_mfields *mflds, int p, int mb, int me);

/* EXTERN_C void cuda_exchange_particles(int p, struct psc_particles *prts); */
/* EXTERN_C void cuda_find_block_indices_ids(struct psc_particles *prts, unsigned int *d_bidx, */
/* 					  unsigned int *d_ids); */
/* EXTERN_C void cuda_find_block_indices_enc_ids(int p, struct psc_particles *prts, unsigned int *d_bidx, */
/* 					      unsigned int *d_ids); */
/* EXTERN_C int  cuda_exclusive_scan(int p, struct psc_particles *prts, unsigned int *d_vals, */
/* 				  unsigned int *d_sums); */
/* EXTERN_C void cuda_reorder_and_offsets(struct psc_particles *prts, unsigned int *d_bidx, unsigned int *d_ids); */

/* EXTERN_C void cuda_sort_patch(int p, struct psc_particles *prts); */
/* EXTERN_C void cuda_sort_patch_by_cell(int p, struct psc_particles *prts); */

/* EXTERN_C void sort_patch_by_cell(int p, struct psc_particles *prts); */

EXTERN_C void sort_pairs_device(unsigned int *d_keys, unsigned int *d_vals, int n);
EXTERN_C void sort_pairs_host(int *d_keys, int *d_vals, int n);
/* EXTERN_C void sort_patch_prep(int p, struct psc_particles *prts, int **d_cnis, */
/* 			      int **d_ids); */
/* EXTERN_C void sort_patch_done(int p, struct psc_particles *prts, int *d_cnis, */
/* 			      int *d_ids); */

/* EXTERN_C void my_sort_check(struct psc_particles *prts, int *keys, int *ids_ref); */

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
