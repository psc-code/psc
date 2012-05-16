
#ifndef PSC_CUDA_H
#define PSC_CUDA_H

#include "psc.h"

#include <assert.h>
#include <math.h>
#include <psc.h>

#include "psc_particles_cuda.h"

// ======================================================================

#define DECLARE_CUDA(pfx)                                               \
  EXTERN_C void pfx##_set_constants(particles_cuda_t *pp,		\
				    fields_cuda_t *pf);			\
  EXTERN_C void pfx##_cuda_push_part_p1(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf,              \
                                        real **d_scratch);              \
  EXTERN_C void pfx##_cuda_push_part_p2(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf);             \
  EXTERN_C void pfx##_cuda_push_part_p3(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf,              \
                                        real *d_scratch,		\
					int block_stride);		\
  EXTERN_C void pfx##_cuda_push_part_p4(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf,              \
                                        real *d_scratch);               \
  EXTERN_C void pfx##_cuda_push_part_p5(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf,              \
                                        real *d_scratch);               \

DECLARE_CUDA(z);
DECLARE_CUDA(z2);
DECLARE_CUDA(z3);
DECLARE_CUDA(yz);
DECLARE_CUDA(yz2);
DECLARE_CUDA(yz3);
DECLARE_CUDA(yz4);
DECLARE_CUDA(yz5);
DECLARE_CUDA(yz6);

EXTERN_C void yz_a_set_constants(particles_cuda_t *pp, fields_cuda_t *pf);
EXTERN_C void yz_b_set_constants(particles_cuda_t *pp, fields_cuda_t *pf);
EXTERN_C void __cuda_push_part_yz_a(particles_cuda_t *pp, fields_cuda_t *pf);
EXTERN_C void __cuda_push_part_yz_b(particles_cuda_t *pp, fields_cuda_t *pf);
EXTERN_C void __cuda_push_part_yz_b3(particles_cuda_t *pp, fields_cuda_t *pf);

struct d_particle {
  real xi[3];
  real qni_div_mni;
  real pxi[3];
  real qni_wni;
};

#define THREADS_PER_BLOCK 256

#define LOAD_PARTICLE(pp, d_p, n) do {					\
    (pp).xi[0]       = d_p.xi4[n].x;					\
    (pp).xi[1]       = d_p.xi4[n].y;					\
    (pp).xi[2]       = d_p.xi4[n].z;					\
    (pp).qni_div_mni = d_p.xi4[n].w;					\
    (pp).pxi[0]      = d_p.pxi4[n].x;					\
    (pp).pxi[1]      = d_p.pxi4[n].y;					\
    (pp).pxi[2]      = d_p.pxi4[n].z;					\
    (pp).qni_wni     = d_p.pxi4[n].w;					\
} while (0)

#define STORE_PARTICLE_POS(pp, d_p, n) do {				\
    d_p.xi4[n].x = (pp).xi[0];						\
    d_p.xi4[n].y = (pp).xi[1];						\
    d_p.xi4[n].z = (pp).xi[2];						\
    d_p.xi4[n].w = (pp).qni_div_mni;					\
} while (0)

#define STORE_PARTICLE_MOM(pp, d_p, n) do {				\
    d_p.pxi4[n].x = (pp).pxi[0];					\
    d_p.pxi4[n].y = (pp).pxi[1];					\
    d_p.pxi4[n].z = (pp).pxi[2];					\
    d_p.pxi4[n].w = (pp).qni_wni;					\
} while (0)

EXTERN_C void __particles_cuda_alloc(particles_cuda_t *pp, bool need_block_offsets,
				     bool need_cell_offsets);
EXTERN_C void __particles_cuda_free(particles_cuda_t *pp);
EXTERN_C void __particles_cuda_to_device(particles_cuda_t *pp,
					 float4 *xi, float4 *pxi,
					 int *offsets, int *c_offsets, int *c_pos);
EXTERN_C void __particles_cuda_to_device_range(particles_cuda_t *pp,
					       float4 *xi, float4 *pxi,
					       int start, int end);
EXTERN_C void __particles_cuda_from_device(particles_cuda_t *pp,
					   float4 *xi4, float4 *pxi4);
EXTERN_C void __particles_cuda_from_device_range(particles_cuda_t *pp,
						 float4 *xi, float4 *pxi,
						 int start, int end);

EXTERN_C void __fields_cuda_alloc(fields_cuda_t *pf);
EXTERN_C void __fields_cuda_free(fields_cuda_t *pf);
EXTERN_C void __fields_cuda_to_device(fields_cuda_t *pf, real *h_flds, int mb, int me);
EXTERN_C void __fields_cuda_from_device(fields_cuda_t *pf, real *h_flds, int mb, int me);

EXTERN_C void cuda_fill_ghosts_periodic_yz(int p, fields_cuda_t *pf, int mb, int me);
EXTERN_C void cuda_fill_ghosts_periodic_z(int p, fields_cuda_t *pf, int mb, int me);
EXTERN_C void cuda_add_ghosts_periodic_yz(int p, fields_cuda_t *pf, int mb, int me);
EXTERN_C void cuda_add_ghosts_periodic_z(int p, fields_cuda_t *pf, int mb, int me);
EXTERN_C void cuda_conducting_wall_E_lo_hi_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_H_lo_hi_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_J_lo_hi_y(int p, fields_cuda_t *pf);

EXTERN_C void cuda_exchange_particles(int p, particles_cuda_t *pp);
EXTERN_C void cuda_alloc_block_indices(particles_cuda_t *pp, unsigned int **d_bidx);
EXTERN_C void cuda_free_block_indices(unsigned int *d_bidx);
EXTERN_C void cuda_find_block_indices_ids(particles_cuda_t *pp, unsigned int *d_bidx,
					  unsigned int *d_ids);
EXTERN_C void cuda_find_block_indices_enc_ids(int p, particles_cuda_t *pp, unsigned int *d_bidx,
					      unsigned int *d_ids);
EXTERN_C int  cuda_exclusive_scan(int p, particles_cuda_t *pp, unsigned int *d_vals,
				  unsigned int *d_sums);
EXTERN_C void cuda_reorder_and_offsets(particles_cuda_t *pp, unsigned int *d_bidx, unsigned int *d_ids);
EXTERN_C void cuda_copy_bidx_from_dev(particles_cuda_t *pp, unsigned int *h_bidx, unsigned int *d_bidx);
EXTERN_C void cuda_copy_offsets_from_dev(particles_cuda_t *pp, unsigned int *h_offsets);

EXTERN_C void cuda_sort_patch(int p, particles_cuda_t *pp);
EXTERN_C void cuda_sort_patch_by_cell(int p, particles_cuda_t *pp);

EXTERN_C void sort_patch_by_cell(int p, particles_cuda_t *pp);

EXTERN_C void sort_pairs_device(unsigned int *d_keys, unsigned int *d_vals, int n);
EXTERN_C void sort_pairs_host(int *d_keys, int *d_vals, int n);
EXTERN_C void sort_patch_prep(int p, particles_cuda_t *pp, int **d_cnis,
			      int **d_ids);
EXTERN_C void sort_patch_done(int p, particles_cuda_t *pp, int *d_cnis,
			      int *d_ids);

EXTERN_C void my_sort_check(particles_cuda_t *pp, int *keys, int *ids_ref);

#define CUDA2_STRIPE_SIZE THREADS_PER_BLOCK

EXTERN_C void psc_mparticles_cuda_get_cuda_2(mparticles_cuda_t *particles,
					     void *particles_base,
					     unsigned int flags);

// These are for field boundary exchange, so they could, eventually, miss the
// interior part

struct cuda_fields_ctx {
  fields_cuda_real_t *arr_off;
  int im[3];
  int ib[3];
  fields_cuda_real_t *arr;
};

struct cuda_mfields_ctx {
  struct cuda_fields_ctx *cf;
};

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
