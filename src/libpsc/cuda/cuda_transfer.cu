
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <mrc_profile.h>

// FIXME, hardcoding is bad, needs to be consistent, etc...
#define BND  (3)
#define MAX_BND_COMPONENTS (3)

EXTERN_C void
__particles_cuda_to_device(struct psc_particles *prts, float4 *xi4, float4 *pxi4)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  struct psc_mparticles *mprts = cuda->mprts;
  assert(mprts);
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  unsigned int off = 0;
  for (int p = 0; p < prts->p; p++) {
    off += psc_mparticles_get_patch(mprts, p)->n_part;
  }
  int n_part = prts->n_part;

  check(cudaMemcpy(mprts_cuda->d_xi4 + off, xi4, n_part * sizeof(*xi4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(mprts_cuda->d_pxi4 + off, pxi4, n_part * sizeof(*pxi4),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__particles_cuda_from_device(struct psc_particles *prts, float4 *xi4, float4 *pxi4)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  struct psc_mparticles *mprts = cuda->mprts;
  assert(mprts);
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  unsigned int off = 0;
  for (int p = 0; p < prts->p; p++) {
    off += psc_mparticles_get_patch(mprts, p)->n_part;
  }
  int n_part = prts->n_part;

  check(cudaMemcpy(xi4, mprts_cuda->d_xi4 + off, n_part * sizeof(*xi4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(pxi4, mprts_cuda->d_pxi4 + off, n_part * sizeof(*pxi4),
		   cudaMemcpyDeviceToHost));
}

EXTERN_C void
cuda_copy_bidx_from_dev(struct psc_particles *prts, unsigned int *h_bidx, unsigned int *d_bidx)
{
  check(cudaMemcpy(h_bidx, d_bidx, prts->n_part * sizeof(*h_bidx),
		   cudaMemcpyDeviceToHost));
}

EXTERN_C void
cuda_copy_bidx_to_dev(struct psc_particles *prts, unsigned int *d_bidx, unsigned int *h_bidx)
{
  check(cudaMemcpy(d_bidx, h_bidx, prts->n_part * sizeof(*d_bidx),
		   cudaMemcpyHostToDevice));
}

void
__psc_mparticles_cuda_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  if (mprts->nr_patches == 0) {
    return;
  }
  
  // FIXME we assume that every patch will have those same dims
  int *ldims = ppsc->patch[0].ldims;

  if (!mprts->flags) {
    // FIXME, they get set too late, so auto-dispatch "1vb" doesn't work
    mprts->flags = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4 | MP_NO_CHECKERBOARD;
  }

  int bs[3];
  for (int d = 0; d < 3; d++) {
    switch (mprts->flags & MP_BLOCKSIZE_MASK) {
    case MP_BLOCKSIZE_1X1X1: bs[d] = 1; break;
    case MP_BLOCKSIZE_2X2X2: bs[d] = 2; break;
    case MP_BLOCKSIZE_4X4X4: bs[d] = 4; break;
    case MP_BLOCKSIZE_8X8X8: bs[d] = 8; break;
    default: assert(0);
    }
    if (ppsc->domain.gdims[d] == 1) {
      bs[d] = 1;
    }
    mprts_cuda->blocksize[d] = bs[d];
    assert(ldims[d] % bs[d] == 0); // FIXME not sure what breaks if not
    mprts_cuda->b_mx[d] = (ldims[d] + bs[d] - 1) / bs[d];
    mprts_cuda->b_dxi[d] = 1.f / (mprts_cuda->blocksize[d] * ppsc->dx[d]);
  }
  mprts_cuda->nr_blocks = mprts_cuda->b_mx[0] * mprts_cuda->b_mx[1] * mprts_cuda->b_mx[2];
  mprts_cuda->nr_total_blocks = mprts->nr_patches * mprts_cuda->nr_blocks;

  mprts_cuda->h_dev = new particles_cuda_dev_t[mprts->nr_patches];
  check(cudaMalloc(&mprts_cuda->d_dev,
		   mprts->nr_patches * sizeof(*mprts_cuda->d_dev)));

  mprts_cuda->nr_prts = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *prts_cuda = psc_particles_cuda(prts);
    mprts_cuda->nr_prts += prts->n_part;
    prts_cuda->mprts = mprts;
  }
  mprts_cuda->h_bnd_cnt = new unsigned int[mprts_cuda->nr_total_blocks];
  unsigned int nr_alloced = mprts_cuda->nr_prts * 1.2;
  mprts_cuda->nr_alloced = nr_alloced;

  check(cudaMalloc((void **) &mprts_cuda->d_xi4, nr_alloced * sizeof(float4)));
  check(cudaMalloc((void **) &mprts_cuda->d_pxi4, nr_alloced * sizeof(float4)));
  check(cudaMalloc((void **) &mprts_cuda->d_alt_xi4, nr_alloced * sizeof(float4)));
  check(cudaMalloc((void **) &mprts_cuda->d_alt_pxi4, nr_alloced * sizeof(float4)));
  check(cudaMalloc((void **) &mprts_cuda->d_bidx, nr_alloced * sizeof(unsigned int)));
  check(cudaMalloc((void **) &mprts_cuda->d_alt_bidx, nr_alloced * sizeof(unsigned int)));
  check(cudaMalloc((void **) &mprts_cuda->d_ids, nr_alloced * sizeof(unsigned int)));
  check(cudaMalloc((void **) &mprts_cuda->d_sums, nr_alloced * sizeof(unsigned int)));

  check(cudaMalloc((void **) &mprts_cuda->d_off, 
		   (mprts_cuda->nr_total_blocks + 1) * sizeof(*mprts_cuda->d_off)));
  check(cudaMalloc((void **) &mprts_cuda->d_bnd_spine_cnts,
		   (1 + mprts_cuda->nr_total_blocks * (CUDA_BND_STRIDE + 1)) * sizeof(unsigned int)));
  check(cudaMalloc((void **) &mprts_cuda->d_bnd_spine_sums,
		   (1 + mprts_cuda->nr_total_blocks * (CUDA_BND_STRIDE + 1)) * sizeof(unsigned int)));

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *prts_cuda = psc_particles_cuda(prts);

    prts_cuda->h_dev = &mprts_cuda->h_dev[p];
    prts_cuda->d_dev = &mprts_cuda->d_dev[p];
  }
}

void
__psc_mparticles_cuda_free(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  delete[] mprts_cuda->h_dev;
  delete[] mprts_cuda->h_bnd_cnt;

  check(cudaFree(mprts_cuda->d_xi4));
  check(cudaFree(mprts_cuda->d_pxi4));
  check(cudaFree(mprts_cuda->d_alt_xi4));
  check(cudaFree(mprts_cuda->d_alt_pxi4));
  check(cudaFree(mprts_cuda->d_bidx));
  check(cudaFree(mprts_cuda->d_alt_bidx));
  check(cudaFree(mprts_cuda->d_ids));
  check(cudaFree(mprts_cuda->d_sums));
  check(cudaFree(mprts_cuda->d_bnd_spine_cnts));
  check(cudaFree(mprts_cuda->d_bnd_spine_sums));

  check(cudaFree(mprts_cuda->d_dev));
}

// ======================================================================
// ======================================================================
// fields

void
__psc_mfields_cuda_setup(struct psc_mfields *mflds)
{
  assert(!ppsc->domain.use_pml);
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  unsigned int total_size = 0;
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
    total_size += size;
  }

  check(cudaMalloc((void **) &mflds_cuda->d_flds,
		   mflds->nr_fields * total_size * sizeof(float)));
  float *d_flds = mflds_cuda->d_flds;

  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);

    unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
    flds_cuda->d_flds = d_flds;
    d_flds += flds->nr_comp * size;
    
    if (flds->im[0] == 1 + 2*BND) {
      int B = 2*BND;
      unsigned int buf_size = 2*B * (flds->im[1] + flds->im[2] - 2*B);
      flds_cuda->h_bnd_buf = new real[MAX_BND_COMPONENTS * buf_size];
      check(cudaMalloc((void **) &flds_cuda->d_bnd_buf,
		       MAX_BND_COMPONENTS * buf_size * sizeof(*flds_cuda->d_bnd_buf)));
    }
  }
}

void
__psc_mfields_cuda_destroy(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  check(cudaFree(mflds_cuda->d_flds));

  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
    
    if (flds->im[0] == 1 + 2*BND) {
      delete[] flds_cuda->h_bnd_buf;
      check(cudaFree(flds_cuda->d_bnd_buf));
    }
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
