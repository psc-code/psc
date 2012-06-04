
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <mrc_profile.h>

EXTERN_C void
cuda_init(int rank)
{
  static bool inited;
  if (!inited) {
    inited = true;
//    cudaSetDevice(rank % 3);
  }
}

// FIXME, hardcoding is bad, needs to be consistent, etc...
#define BND  (3)
#define MAX_BND_COMPONENTS (3)

EXTERN_C void
__particles_cuda_alloc(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  int n_alloced = prts->n_part * 1.2; // FIXME, need to handle realloc eventualy
  cuda->n_alloced = n_alloced;
  particles_cuda_dev_t *h_dev = cuda->h_dev;

  check(cudaMalloc((void **) &h_dev->xi4, n_alloced * sizeof(float4)));
  check(cudaMalloc((void **) &h_dev->pxi4, n_alloced * sizeof(float4)));

  check(cudaMalloc((void **) &h_dev->alt_xi4, n_alloced * sizeof(float4)));
  check(cudaMalloc((void **) &h_dev->alt_pxi4, n_alloced * sizeof(float4)));

  check(cudaMalloc((void **) &h_dev->offsets, 
		   (cuda->nr_blocks + 1) * sizeof(int)));
  check(cudaMemcpy(&h_dev->offsets[cuda->nr_blocks], &prts->n_part, sizeof(int),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__particles_cuda_to_device(struct psc_particles *prts, float4 *xi4, float4 *pxi4,
			   int *offsets)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  int n_part = prts->n_part;
  particles_cuda_dev_t *h_dev = cuda->h_dev;

  assert(n_part <= cuda->n_alloced);
  check(cudaMemcpy(h_dev->xi4, xi4, n_part * sizeof(*xi4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(h_dev->pxi4, pxi4, n_part * sizeof(*pxi4),
		   cudaMemcpyHostToDevice));
  if (offsets) {
    check(cudaMemcpy(h_dev->offsets, offsets,
                    (cuda->nr_blocks + 1) * sizeof(int), cudaMemcpyHostToDevice));
  }
}

EXTERN_C void
__particles_cuda_to_device_range(struct psc_particles *prts, float4 *xi4, float4 *pxi4,
				 int start, int end)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  particles_cuda_dev_t *h_dev = cuda->h_dev;

  check(cudaMemcpy(h_dev->xi4 + start, xi4, (end - start) * sizeof(*xi4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(h_dev->pxi4 + start, pxi4, (end - start) * sizeof(*pxi4),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__particles_cuda_from_device(struct psc_particles *prts, float4 *xi4, float4 *pxi4)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  int n_part = prts->n_part;
  particles_cuda_dev_t *h_dev = cuda->h_dev;

  assert(n_part <= cuda->n_alloced);
  check(cudaMemcpy(xi4, h_dev->xi4, n_part * sizeof(*xi4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(pxi4, h_dev->pxi4, n_part * sizeof(*pxi4),
		   cudaMemcpyDeviceToHost));
}

EXTERN_C void
__particles_cuda_from_device_range(struct psc_particles *prts, float4 *xi4, float4 *pxi4,
				   int start, int end)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  particles_cuda_dev_t *h_dev = cuda->h_dev;

  check(cudaMemcpy(xi4, h_dev->xi4 + start, (end - start) * sizeof(*xi4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(pxi4, h_dev->pxi4 + start, (end - start) * sizeof(*pxi4),
		   cudaMemcpyDeviceToHost));
}

EXTERN_C void
__particles_cuda_free(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  particles_cuda_dev_t *h_dev = cuda->h_dev;

  check(cudaFree(h_dev->xi4));
  check(cudaFree(h_dev->pxi4));
  check(cudaFree(h_dev->alt_xi4));
  check(cudaFree(h_dev->alt_pxi4));
  check(cudaFree(h_dev->offsets));
}

EXTERN_C void
cuda_copy_offsets_from_dev(struct psc_particles *prts, unsigned int *h_offsets)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  check(cudaMemcpy(h_offsets, cuda->h_dev->offsets, (cuda->nr_blocks + 1) * sizeof(*h_offsets),
		   cudaMemcpyDeviceToHost));
}

EXTERN_C void
cuda_copy_offsets_to_dev(struct psc_particles *prts, unsigned int *h_offsets)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  check(cudaMemcpy(cuda->h_dev->offsets, h_offsets, (cuda->nr_blocks + 1) * sizeof(*h_offsets),
		   cudaMemcpyHostToDevice));
}

void
__psc_mparticles_cuda_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  mprts_cuda->h_dev = new particles_cuda_dev_t[mprts->nr_patches];
  check(cudaMalloc(&mprts_cuda->d_dev,
		   mprts->nr_patches * sizeof(*mprts_cuda->d_dev)));

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *prts_cuda = psc_particles_cuda(prts);
    prts_cuda->h_dev = &mprts_cuda->h_dev[p];

    __particles_cuda_alloc(prts);
    cuda_alloc_block_indices(prts, &prts_cuda->h_dev->bidx); // FIXME, merge into ^^^
    cuda_alloc_block_indices(prts, &prts_cuda->h_dev->ids);
    cuda_alloc_block_indices(prts, &prts_cuda->h_dev->alt_bidx);
    cuda_alloc_block_indices(prts, &prts_cuda->h_dev->alt_ids);
    cuda_alloc_block_indices(prts, &prts_cuda->h_dev->sums);

    prts_cuda->d_dev = &mprts_cuda->d_dev[p];
  }
}

void
__psc_mparticles_cuda_free(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *prts_cuda = psc_particles_cuda(prts);
    cuda_free_block_indices(prts_cuda->h_dev->bidx);
    cuda_free_block_indices(prts_cuda->h_dev->ids);
    cuda_free_block_indices(prts_cuda->h_dev->alt_bidx);
    cuda_free_block_indices(prts_cuda->h_dev->alt_ids);
    cuda_free_block_indices(prts_cuda->h_dev->sums);
    __particles_cuda_free(prts);
  }
  free(mprts_cuda->h_dev);
  check(cudaFree(mprts_cuda->d_dev));
}

// ======================================================================
// ======================================================================
// fields

EXTERN_C void
__fields_cuda_alloc(struct psc_fields *pf)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  assert(!ppsc->domain.use_pml);

  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMalloc((void **) &pfc->d_flds, pf->nr_comp * size * sizeof(float)));

  if (pf->im[0] == 1 + 2*BND) {
    int B = 2*BND;
    unsigned int buf_size = 2*B * (pf->im[1] + pf->im[2] - 2*B);
    pfc->h_bnd_buf = new real[MAX_BND_COMPONENTS * buf_size];
    check(cudaMalloc((void **) &pfc->d_bnd_buf,
		     MAX_BND_COMPONENTS * buf_size * sizeof(*pfc->d_bnd_buf)));
  }
}

EXTERN_C void
__fields_cuda_free(struct psc_fields *pf)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  check(cudaFree(pfc->d_flds));

  if (pf->im[0] == 1 + 2*BND) {
    delete[] pfc->h_bnd_buf;
    check(cudaFree(pfc->d_bnd_buf));
  }
}

EXTERN_C void
__fields_cuda_to_device(struct psc_fields *pf, real *h_flds, int mb, int me)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMemcpy(pfc->d_flds + mb * size,
		   h_flds + mb * size,
		   (me - mb) * size * sizeof(float),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__fields_cuda_from_device(struct psc_fields *pf, real *h_flds, int mb, int me)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMemcpy(h_flds + mb * size,
		   pfc->d_flds + mb * size,
		   (me - mb) * size * sizeof(float),
		   cudaMemcpyDeviceToHost));
}

// ======================================================================

enum {
  PACK,
  UNPACK,
};

// ======================================================================
// fields_device_pack

// FIXME/OPT: can probably be accelerated by making component the fast index

template<int B, int what>
__global__ static void
k_fields_device_pack_yz(real *d_buf, real *d_flds, int gmy, int gmz, int mm)
{
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);
  int gmx = 2*BND + 1;
  int jx = BND;
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  int n_threads = mm * buf_size;
  if (tid >= n_threads)
    return;

  int n = tid;
  int m = n / buf_size; n -= m * buf_size;
  int jz, jy;
  if (n < B * gmy) {
    jz = n / gmy; n -= jz * gmy;
    jy = n;
  } else if (n < B * gmy + (gmz - 2*B) * 2*B) {
    n -= B * gmy;
    jz = n / (2*B); n -= jz * 2*B;
    if (n < B) {
      jy = n;
    } else {
      jy = n + gmy - 2*B;
    }
    jz += B;
  } else {
    n -= B * gmy + (gmz - 2*B) * 2*B;
    jz = n / gmy; n -= jz * gmy;
    jy = n;
    jz += gmz - B;
  }
  
  // FIXME, should use F3_DEV_YZ
  if (what == PACK) {
    d_buf[tid] = d_flds[((m * gmz + jz) * gmy + jy) * gmx + jx];
  } else if (what == UNPACK) {
    d_flds[((m * gmz + jz) * gmy + jy) * gmx + jx] = d_buf[tid]; 
  }
}

template<int B, bool pack>
static void
fields_device_pack_yz(struct psc_fields *pf, int mb, int me)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  int gmy = pf->im[1], gmz = pf->im[2];
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);
  int n_threads = buf_size * (me - mb);

  dim3 dimGrid((n_threads + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_fields_device_pack_yz<B, pack> <<<dimGrid, dimBlock>>>
    (pfc->d_bnd_buf, pfc->d_flds + mb * size, gmy, gmz, me - mb);
}

// ======================================================================
// fields_host_pack

#define WHAT do {					\
    if (what == PACK) {					\
      h_buf[tid++] = F3_CF_0(cf, m, 0,jy,jz);		\
    } else if (what == UNPACK) {			\
      F3_CF_0(cf, m, 0,jy,jz) = h_buf[tid++];		\
    }							\
  } while(0)

template<int B, int what>
static void
fields_host_pack_yz(struct psc_fields_cuda_bnd *cf, real *h_buf, int mb, int me)
{
  int gmy = cf->im[1], gmz = cf->im[2];
  int tid = 0;
  for (int m = 0; m < me - mb; m++) {
    for (int jz = 0; jz < B; jz++) {
      for (int jy = 0; jy < gmy; jy++) {
	WHAT;
      }
    }
    for (int jz = B; jz < gmz - B; jz++) {
      for (int jy = 0; jy < B; jy++) {
	WHAT;
      }
      for (int jy = gmy - B; jy < gmy; jy++) {
	WHAT;
      }
    }
    for (int jz = gmz - B; jz < gmz; jz++) {
      for (int jy = 0; jy < gmy; jy++) {
	WHAT;
      }
    }
  }
}

template<int B>
static void
__fields_cuda_from_device_yz(struct psc_fields *pf, struct psc_fields_cuda_bnd *cf, int mb, int me)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  int gmy = pf->im[1], gmz = pf->im[2];
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);
  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(pf->ib[1] == -BND);
  assert(pf->im[1] >= 2 * B);
  assert(pf->im[2] >= 2 * B);

  // static int pr1, pr2, pr3;
  // if (!pr1) {
  //   pr1 = prof_register("field_device_pack", 1., 0, 0);
  //   pr2 = prof_register("cuda_memcpy", 1., 0, 0);
  //   pr3 = prof_register("field_host_unpack", 1., 0, 0);
  // }

//  prof_start(pr1);
  fields_device_pack_yz<B, PACK>(pf, mb, me);
  cuda_sync_if_enabled();
//  prof_stop(pr1);
//  prof_start(pr2);
  check(cudaMemcpy(pfc->h_bnd_buf, pfc->d_bnd_buf,
		   (me - mb) * buf_size * sizeof(*pfc->h_bnd_buf), cudaMemcpyDeviceToHost));
//  prof_stop(pr2);
//  prof_start(pr3);
  fields_host_pack_yz<B, UNPACK>(cf, pfc->h_bnd_buf, mb, me);
//  prof_stop(pr3);
}

template<int B>
static void
__fields_cuda_to_device_yz(struct psc_fields *pf, struct psc_fields_cuda_bnd *cf, int mb, int me)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  int gmy = pf->im[1], gmz = pf->im[2];
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);
  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(pf->ib[1] == -BND);
  assert(pf->im[1] >= 2 * B);
  assert(pf->im[2] >= 2 * B);

  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_host_pack", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy", 1., 0, 0);
    pr3 = prof_register("field_device_unpack", 1., 0, 0);
  }

  prof_start(pr1);
  fields_host_pack_yz<B, PACK>(cf, pfc->h_bnd_buf, mb, me);
  prof_stop(pr1);
  prof_start(pr2);
  check(cudaMemcpy(pfc->d_bnd_buf, pfc->h_bnd_buf,
		   (me - mb) * buf_size * sizeof(*pfc->d_bnd_buf), cudaMemcpyHostToDevice));
  prof_stop(pr2);
  prof_start(pr3);
  fields_device_pack_yz<B, UNPACK>(pf, mb, me);
  cuda_sync_if_enabled();
  prof_stop(pr3);
}

// ======================================================================

EXTERN_C void
__fields_cuda_from_device_inside(struct psc_fields *pf, int mb, int me)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  if (pf->im[0] == 2 * -pf->ib[0] + 1) {
    __fields_cuda_from_device_yz<2*BND>(pf, &pfc->bnd, mb, me);
  } else {
    unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
    check(cudaMemcpy(pfc->bnd.arr,
		     pfc->d_flds + mb * size,
		     (me - mb) * size * sizeof(float),
		     cudaMemcpyDeviceToHost));
  }
}

EXTERN_C void
__fields_cuda_to_device_outside(struct psc_fields *pf, int mb, int me)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  if (pf->im[0] == 2 * -pf->ib[0] + 1) {
    __fields_cuda_to_device_yz<BND>(pf, &pfc->bnd, mb, me);
  } else {
    unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
    check(cudaMemcpy(pfc->d_flds + mb * size,
		     pfc->bnd.arr,
		     (me - mb) * size * sizeof(float),
		     cudaMemcpyHostToDevice));
  }
}

EXTERN_C void
__fields_cuda_to_device_inside(struct psc_fields *pf, int mb, int me)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  if (pf->im[0] == 2 * -pf->ib[0] + 1) {
    __fields_cuda_to_device_yz<2*BND>(pf, &pfc->bnd, mb, me);
  } else {
    unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
    check(cudaMemcpy(pfc->d_flds + mb * size,
		     pfc->bnd.arr,
		     (me - mb) * size * sizeof(float),
		     cudaMemcpyHostToDevice));
  }
}
