
#include "cuda_mparticles.h"

#include "psc_cuda.h"
#include "particles_cuda.h"

#include <mrc_ddc_private.h>
#include <mrc_profile.h>

void
psc_cuda_init(void)
{
  static bool first_time = true;
  if (!first_time)
    return;

  first_time = false;

  int deviceCount;
  cudaGetDeviceCount(&deviceCount);

  // This function call returns 0 if there are no CUDA capable devices.
  if (deviceCount == 0) {
    printf("There is no device supporting CUDA\n");
    return;
  }

  for (int dev = 0; dev < deviceCount; ++dev) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);

    if (dev == 0) {
      // This function call returns 9999 for both major & minor fields, if no CUDA capable devices are present
      if (deviceProp.major == 9999 && deviceProp.minor == 9999)
	printf("There is no device supporting CUDA.\n");
      else if (deviceCount == 1)
	printf("There is 1 device supporting CUDA\n");
      else
	printf("There are %d devices supporting CUDA\n", deviceCount);
    }
    printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
    printf("  CUDA Capability Major revision number:         %d\n", deviceProp.major);
    printf("  CUDA Capability Minor revision number:         %d\n", deviceProp.minor);
    printf("  Total amount of global memory:                 %lu bytes\n", deviceProp.totalGlobalMem);
#if CUDART_VERSION >= 2000
    printf("  Number of multiprocessors:                     %d\n", deviceProp.multiProcessorCount);
    printf("  Number of cores:                               %d\n", 8 * deviceProp.multiProcessorCount);
#endif
    printf("  Total amount of constant memory:               %lu bytes\n", deviceProp.totalConstMem); 
    printf("  Total amount of shared memory per block:       %lu bytes\n", deviceProp.sharedMemPerBlock);
    printf("  Total number of registers available per block: %d\n", deviceProp.regsPerBlock);
    printf("  Warp size:                                     %d\n", deviceProp.warpSize);
    printf("  Maximum number of threads per block:           %d\n", deviceProp.maxThreadsPerBlock);
    printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
	   deviceProp.maxThreadsDim[0],
	   deviceProp.maxThreadsDim[1],
	   deviceProp.maxThreadsDim[2]);
    printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
	   deviceProp.maxGridSize[0],
	   deviceProp.maxGridSize[1],
	   deviceProp.maxGridSize[2]);
    printf("  Maximum memory pitch:                          %lu bytes\n", deviceProp.memPitch);
    printf("  Texture alignment:                             %lu bytes\n", deviceProp.textureAlignment);
    printf("  Clock rate:                                    %.2f GHz\n", deviceProp.clockRate * 1e-6f);
#if CUDART_VERSION >= 2000
    printf("  Concurrent copy and execution:                 %s\n", deviceProp.deviceOverlap ? "Yes" : "No");
#endif
#if CUDART_VERSION >= 2020
    printf("  Run time limit on kernels:                     %s\n", deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
    printf("  Integrated:                                    %s\n", deviceProp.integrated ? "Yes" : "No");
    printf("  Support host page-locked memory mapping:       %s\n", deviceProp.canMapHostMemory ? "Yes" : "No");
    printf("  Compute mode:                                  %s\n", deviceProp.computeMode == cudaComputeModeDefault ?
	   "Default (multiple host threads can use this device simultaneously)" :
	   deviceProp.computeMode == cudaComputeModeExclusive ?
	   "Exclusive (only one host thread at a time can use this device)" :
	   deviceProp.computeMode == cudaComputeModeProhibited ?
	   "Prohibited (no host thread can use this device)" :
	   "Unknown");
#endif
  }
}
	   
	   
// FIXME, hardcoding is bad, needs to be consistent, etc...
#define MAX_BND_COMPONENTS (3)

EXTERN_C void
__particles_cuda_to_device(struct psc_particles *prts, float4 *xi4, float4 *pxi4)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  struct psc_mparticles *mprts = cuda->mprts;
  assert(mprts);
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;
  assert(cmprts);
  unsigned int off = 0;
  for (int p = 0; p < prts->p; p++) {
    off += psc_mparticles_get_patch(mprts, p)->n_part;
  }
  int n_part = prts->n_part;

  check(cudaMemcpy(cmprts->d_xi4 + off, xi4, n_part * sizeof(*xi4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(cmprts->d_pxi4 + off, pxi4, n_part * sizeof(*pxi4),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__particles_cuda_from_device(struct psc_particles *prts, float4 *xi4, float4 *pxi4)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  struct psc_mparticles *mprts = cuda->mprts;
  assert(mprts);
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;
  assert(cmprts);
  unsigned int off = 0;
  for (int p = 0; p < prts->p; p++) {
    off += psc_mparticles_get_patch(mprts, p)->n_part;
  }
  int n_part = prts->n_part;

  check(cudaMemcpy(xi4, cmprts->d_xi4 + off, n_part * sizeof(*xi4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(pxi4, cmprts->d_pxi4 + off, n_part * sizeof(*pxi4),
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
  psc_cuda_init();

  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = cuda_mparticles_create();
  mprts_cuda->cmprts = cmprts;

  cmprts->n_patches = mprts->nr_patches;

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
    // assumes no AMR
    mprts_cuda->b_dxi[d] = 1.f / (mprts_cuda->blocksize[d] * ppsc->patch[0].dx[d]);
  }
  mprts_cuda->nr_blocks = mprts_cuda->b_mx[0] * mprts_cuda->b_mx[1] * mprts_cuda->b_mx[2];
  mprts_cuda->nr_total_blocks = mprts->nr_patches * mprts_cuda->nr_blocks;

  mprts_cuda->h_n_prts = new int[mprts->nr_patches];
  check(cudaMalloc(&mprts_cuda->d_n_prts,
		   mprts->nr_patches * sizeof(int)));

  cmprts->n_prts = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *prts_cuda = psc_particles_cuda(prts);
    cmprts->n_prts += prts->n_part;
    prts_cuda->mprts = mprts;
  }
  mprts_cuda->h_bnd_cnt = new unsigned int[mprts_cuda->nr_total_blocks];
  unsigned int nr_alloced = cmprts->n_prts * 1.4;
  mprts_cuda->nr_alloced = nr_alloced;

  cuda_mparticles_alloc(cmprts, nr_alloced);

  check(cudaMalloc((void **) &mprts_cuda->d_alt_bidx, nr_alloced * sizeof(unsigned int)));
  check(cudaMalloc((void **) &mprts_cuda->d_sums, nr_alloced * sizeof(unsigned int)));

  check(cudaMalloc((void **) &mprts_cuda->d_off, 
		   (mprts_cuda->nr_total_blocks + 1) * sizeof(*mprts_cuda->d_off)));
  check(cudaMalloc((void **) &mprts_cuda->d_bnd_spine_cnts,
		   (1 + mprts_cuda->nr_total_blocks * (CUDA_BND_STRIDE + 1)) * sizeof(unsigned int)));
  check(cudaMalloc((void **) &mprts_cuda->d_bnd_spine_sums,
		   (1 + mprts_cuda->nr_total_blocks * (CUDA_BND_STRIDE + 1)) * sizeof(unsigned int)));
}

void
__psc_mparticles_cuda_free(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;
  assert(cmprts);

  delete[] mprts_cuda->h_n_prts;
  delete[] mprts_cuda->h_bnd_cnt;

  cuda_mparticles_dealloc(cmprts);

  check(cudaFree(mprts_cuda->d_alt_bidx));
  check(cudaFree(mprts_cuda->d_sums));
  check(cudaFree(mprts_cuda->d_off));
  check(cudaFree(mprts_cuda->d_bnd_spine_cnts));
  check(cudaFree(mprts_cuda->d_bnd_spine_sums));

  check(cudaFree(mprts_cuda->d_n_prts));

  cuda_mparticles_destroy(cmprts);
  mprts_cuda->cmprts = NULL;
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
  unsigned int buf_size = 0;
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    if (p == 0) {
      for (int d = 0; d < 3; d++) {
	mflds_cuda->im[d] = flds->im[d];
	mflds_cuda->ib[d] = flds->ib[d];
      }
    } else {
      for (int d = 0; d < 3; d++) {
	assert(mflds_cuda->im[d] == flds->im[d]);
	assert(mflds_cuda->ib[d] == flds->ib[d]);
      }
    }
    unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
    total_size += size;
    if (flds->im[0] == 1) {;// + 2*BND) {
      int B = 2*BND;
      buf_size = 2*B * (flds->im[1] + flds->im[2] - 2*B);
    } else {
      assert(0);
    }
  }

  mprintf("nr_fields %d tsize %d\n", mflds->nr_fields, total_size);
  check(cudaMalloc((void **) &mflds_cuda->d_flds,
		   mflds->nr_fields * total_size * sizeof(float)));
  check(cudaMalloc((void **) &mflds_cuda->d_bnd_buf,
		   MAX_BND_COMPONENTS * buf_size * mflds->nr_patches * sizeof(float)));
  mflds_cuda->h_bnd_buf = new float[MAX_BND_COMPONENTS * mflds->nr_patches * buf_size];
  float *d_flds = mflds_cuda->d_flds;

  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    assert(psc_fields_ops(flds) == &psc_fields_cuda_ops);
    struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);

    unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
    flds_cuda->d_flds = d_flds;
    assert(d_flds == mflds_cuda->d_flds + p * flds->nr_comp * size);
    d_flds += flds->nr_comp * size;
    
    struct psc_fields_cuda_bnd *cf = &flds_cuda->bnd;
    int sz = 1;
    for (int d = 0; d < 3; d++) {
      if (flds->im[d] == 1 - 2 * flds->ib[d]) { // only 1 non-ghost point
	cf->im[d] = 1;
	cf->ib[d] = 0;
      } else {
	cf->im[d] = flds->im[d];
	cf->ib[d] = flds->ib[d];
      }
      sz *= cf->im[d];
    }
    cf->arr = new float [MAX_BND_COMPONENTS * sz];
    cf->arr_off = cf->arr 
      - ((cf->ib[2] * cf->im[1] + cf->ib[1]) * cf->im[0] + cf->ib[0]);
  }
}

void
__psc_mfields_cuda_destroy(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  check(cudaFree(mflds_cuda->d_flds));
  check(cudaFree(mflds_cuda->d_bnd_buf));
  check(cudaFree(mflds_cuda->d_nei_patch));
  check(cudaFree(mflds_cuda->d_map_out));
  check(cudaFree(mflds_cuda->d_map_in));
  delete[] mflds_cuda->h_bnd_buf;
  delete[] mflds_cuda->h_nei_patch;
  delete[] mflds_cuda->h_map_out;
  delete[] mflds_cuda->h_map_in;

  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
    struct psc_fields_cuda_bnd *cf = &flds_cuda->bnd;
    delete[] cf->arr;
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

template<int B, int WHAT, int NR_COMPONENTS>
__global__ static void
k_fields_device_pack_yz(real *d_buf, real *d_flds, int gmy, int gmz,
			int nr_patches, int nr_fields)
{
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);
  int gmx = 1;//2*BND + 1;
  int jx = 0;//BND;
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  int n_threads = NR_COMPONENTS * buf_size;
  int p = tid / n_threads;
  if (p >= nr_patches)
    return;

  int n = tid - p * n_threads;
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
  if (WHAT == PACK) {
    d_buf[tid] = d_flds[(((p * nr_fields + m) * gmz + jz) * gmy + jy) * gmx + jx];
  } else if (WHAT == UNPACK) {
    d_flds[(((p * nr_fields + m) * gmz + jz) * gmy + jy) * gmx + jx] = d_buf[tid]; 
  }
}

template<int B, int WHAT, int NR_COMPONENTS>
__global__ static void
k_fields_device_pack2_yz(real *d_buf, real *d_flds, int *d_nei_patch_by_dir1,
			 int gmy, int gmz, int nr_patches, int nr_fields)
{
  unsigned int nr_ghosts = 2*B * (gmy + gmz - 2*B);
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  int p = tid / nr_ghosts;
  if (p >= nr_patches)
    return;

  int n = tid - p * nr_ghosts;
  int jy, jz;
  int diry, dirz;
  if (n < 2*B * gmy) {
    jz = n / gmy;
    jy = n - jz * gmy;
    if (jy < B) {
      diry = -1;
    } else if (jy < gmy - B) {
      diry = 0;
    } else {
      diry = 1;
    }
    if (jz < B) {
      dirz = -1;
    } else {
      dirz = 1;
      jz += gmz - 2*B;
    }
  } else {
    n -= 2*B * gmy;
    jz = n / (2*B) + B;
    jy = n % (2*B);
    dirz = 0;
    if (jy < B) {
      diry = -1;
    } else {
      diry = 1;
      jy += gmy - 2*B;
    }
  }

  int s_p = d_nei_patch_by_dir1[p*9 + 3*dirz + diry + 4];
  // copy only ghost areas that interface with remote patches
  if (1||s_p < 0) {
    for (int m = 0; m < NR_COMPONENTS; m++) {
      // FIXME, should use F3_DEV_YZ
      if (WHAT == PACK) {
	d_buf[m * nr_ghosts * nr_patches + tid] = d_flds[((p * nr_fields + m) * gmz + jz) * gmy + jy];
      } else if (WHAT == UNPACK) {
	d_flds[((p * nr_fields + m) * gmz + jz) * gmy + jy] = d_buf[m * nr_ghosts * nr_patches + tid]; 
      }
    }
  }
}

template<int B, int NR_COMPONENTS>
__global__ static void
k_fill_ghosts_local_yz(float *d_flds, int *d_nei_patch_by_dir1,
		       int gmy, int gmz, int nr_fields, int nr_patches)
{
  int nr_ghosts = 2*B * (gmy + gmz - 2*B);
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  int r_p = tid / nr_ghosts;
  if (r_p >= nr_patches)
    return;

  int n = tid - r_p * nr_ghosts;

  int jy, jz;
  int diry, dirz;
  if (n < 2*B * gmy) {
    jz = n / gmy;
    jy = n - jz * gmy;
    if (jy < B) {
      diry = -1;
    } else if (jy < gmy - B) {
      diry = 0;
    } else {
      diry = 1;
    }
    if (jz < B) {
      dirz = -1;
    } else {
      dirz = 1;
      jz += gmz - 2*B;
    }
  } else {
    n -= 2*B * gmy;
    jz = n / (2*B) + B;
    jy = n % (2*B);
    dirz = 0;
    if (jy < B) {
      diry = -1;
    } else {
      diry = 1;
      jy += gmy - 2*B;
    }
  }
  int s_p = d_nei_patch_by_dir1[r_p*9 + 3*dirz + diry + 4];
  if (s_p >= 0) {
    float *r_f = &d_flds[((r_p * nr_fields) * gmz) * gmy];
    float *s_f = &d_flds[((s_p * nr_fields)
			  * gmz - dirz * (gmz - 2*2)) 
			 * gmy - diry * (gmy - 2*2)];
    for (int m = 0; m < NR_COMPONENTS; m++) {
      int i = (m * gmz + jz) * gmy + jy;
      r_f[i] = s_f[i];
    }
  }
}

template<int B, bool pack>
static void
fields_device_pack_yz(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  unsigned int size = mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];
  int gmy = mflds_cuda->im[1], gmz = mflds_cuda->im[2];
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);
  int n_threads = buf_size * (me - mb) * mflds->nr_patches;

  dim3 dimGrid((n_threads + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  float *d_bnd_buf = mflds_cuda->d_bnd_buf;
  float *d_flds = mflds_cuda->d_flds + mb * size;
  if (me - mb == 3) {
    k_fields_device_pack_yz<B, pack, 3> <<<dimGrid, dimBlock>>>
      (d_bnd_buf, d_flds, gmy, gmz, mflds->nr_patches,
       mflds->nr_fields);
  } else if (me - mb == 1) {
    k_fields_device_pack_yz<B, pack, 1> <<<dimGrid, dimBlock>>>
      (d_bnd_buf, d_flds, gmy, gmz, mflds->nr_patches,
       mflds->nr_fields);
  } else {
    assert(0);
  }
  cuda_sync_if_enabled();
}

template<int B, bool pack>
static void
fields_device_pack2_yz(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  const int NR_COMPONENTS = 3;
  assert(me - mb == NR_COMPONENTS);

  int *im = mflds_cuda->im;
  assert(im[0] == 1);
  int nr_fields = mflds->nr_fields;
  int nr_patches = mflds->nr_patches;
  int nr_ghosts = 2*B * (im[1] + im[2] - 2*B);
  int nr_threads = nr_ghosts * nr_patches;

  dim3 dimGrid((nr_threads + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  float *d_flds = mflds_cuda->d_flds + mb * im[1] * im[2];
  k_fields_device_pack2_yz<B, pack, NR_COMPONENTS> <<<dimGrid, dimBlock>>>
    (mflds_cuda->d_bnd_buf, d_flds, mflds_cuda->d_nei_patch,
     im[1], im[2], nr_patches, nr_fields);
  cuda_sync_if_enabled();
}

EXTERN_C void
__fields_cuda_fill_ghosts_local(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  const int B = 2;

  int *im = mflds_cuda->im;
  assert(im[0] == 1);
  int nr_fields = mflds->nr_fields;
  int nr_patches = mflds->nr_patches;
  int nr_ghosts = 2*B * (im[1] + im[2] - 2*B);
  int nr_threads = nr_ghosts * nr_patches;

  dim3 dimGrid((nr_threads + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  float *d_flds = mflds_cuda->d_flds + mb * im[1] * im[2];
  if (me - mb == 3) {
    k_fill_ghosts_local_yz<B, 3> <<<dimGrid, dimBlock>>>
      (d_flds, mflds_cuda->d_nei_patch, im[1], im[2],
       nr_fields, nr_patches);
  } else if (me - mb == 1) {
    k_fill_ghosts_local_yz<B, 1> <<<dimGrid, dimBlock>>>
      (d_flds, mflds_cuda->d_nei_patch, im[1], im[2],
       nr_fields, nr_patches);
  } else {
    assert(0);
  }
  cuda_sync_if_enabled();

#if 0
  thrust::device_ptr<float> d_flds(mflds_cuda->d_flds);
  thrust::host_vector<float> h_flds(d_flds, d_flds + nr_patches * nr_fields * im[2] *im[2]);

  for (int tid = 0; tid < nr_threads; tid++) {
    cuda_fill_ghosts_local_gold(&h_flds[0], nei_patch_by_dir1, mb, me, im, nr_fields, 
				nr_patches, nr_ghosts, tid);
  }
  
  thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
#endif
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
fields_host_pack_yz(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  int gmy = mflds_cuda->im[1], gmz = mflds_cuda->im[2];
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);

  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
    struct psc_fields_cuda_bnd *cf = &flds_cuda->bnd;
    real *h_buf = mflds_cuda->h_bnd_buf + p * buf_size * (me - mb);
    
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
}

template<int B, int what>
static void
fields_host_pack2_yz(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  real *h_buf = mflds_cuda->h_bnd_buf;
  int tid = 0;
  for (int m = 0; m < me - mb; m++) {
    for (int p = 0; p < mflds->nr_patches; p++) {
      struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
      struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
      struct psc_fields_cuda_bnd *cf = &flds_cuda->bnd;
      
      int gmy = cf->im[1], gmz = cf->im[2];
      for (int jz = 0; jz < B; jz++) {
	for (int jy = 0; jy < gmy; jy++) {
	  WHAT;
	}
      }
      for (int jz = gmz - B; jz < gmz; jz++) {
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
    }
  }
}

#undef WHAT

template<int B, int what>
static void
fields_host_pack3_yz(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  int *im = mflds_cuda->im;
  int nr_fields = mflds->nr_fields;

  real *h_buf = mflds_cuda->h_bnd_buf;
  int nr_map;
  int *h_map;
  if (B == 2) {
    h_map = mflds_cuda->h_map_out;
    nr_map = mflds_cuda->nr_map_out;
  } else if (B == 4) {
    h_map = mflds_cuda->h_map_in;
    nr_map = mflds_cuda->nr_map_in;
  } else {
    assert(0);
  }
  for (int tid = 0; tid < nr_map; tid++) {
    int i = h_map[tid];
    int p = i / (nr_fields * im[2] * im[1]);
    int off = i - p * (nr_fields * im[2] * im[1]);
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
    struct psc_fields_cuda_bnd *cf = &flds_cuda->bnd;
    for (int m = 0; m < me - mb; m++) {
      if (what == PACK) {
	h_buf[tid + m * nr_map] = cf->arr[off + m * im[2] * im[1]];
      } else if (what == UNPACK) {
	cf->arr[off + m * im[2] * im[1]] = h_buf[tid + m * nr_map];
      }
    }
  }
}

template<int B>
static void
__fields_cuda_from_device_yz(struct psc_mfields *mflds, int mb, int me)
{
  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_device_pack", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy", 1., 0, 0);
    pr3 = prof_register("field_host_unpack", 1., 0, 0);
  }
    
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  int gmy = mflds_cuda->im[1], gmz = mflds_cuda->im[2];
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);
  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(mflds_cuda->ib[1] == -BND);
  assert(mflds_cuda->im[1] >= 2 * B);
  assert(mflds_cuda->im[2] >= 2 * B);

  prof_start(pr1);
  fields_device_pack_yz<B, PACK>(mflds, mb, me);
  prof_stop(pr1);
  
  prof_start(pr2);
  check(cudaMemcpy(mflds_cuda->h_bnd_buf, mflds_cuda->d_bnd_buf,
		   (me - mb) * buf_size * mflds->nr_patches *
		   sizeof(*mflds_cuda->h_bnd_buf),
		   cudaMemcpyDeviceToHost));
  prof_stop(pr2);

  prof_start(pr3);
  fields_host_pack_yz<B, UNPACK>(mflds, mb, me);
  prof_stop(pr3);
}

template<int B>
static void
__fields_cuda_to_device_yz(struct psc_mfields *mflds, int mb, int me)
{
  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_host_pack", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy", 1., 0, 0);
    pr3 = prof_register("field_device_unpack", 1., 0, 0);
  }
  
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  int gmy = mflds_cuda->im[1], gmz = mflds_cuda->im[2];
  unsigned int buf_size = 2*B * (gmy + gmz - 2*B);
  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(mflds_cuda->ib[1] == -BND);
  assert(mflds_cuda->im[1] >= 2 * B);
  assert(mflds_cuda->im[2] >= 2 * B);

  prof_start(pr1);
  fields_host_pack_yz<B, PACK>(mflds, mb, me);
  prof_stop(pr1);

  prof_start(pr2);
  check(cudaMemcpy(mflds_cuda->d_bnd_buf, mflds_cuda->h_bnd_buf,
		   MAX_BND_COMPONENTS * buf_size * mflds->nr_patches *
		   sizeof(*mflds_cuda->d_bnd_buf),
		   cudaMemcpyHostToDevice));
  prof_stop(pr2);

  prof_start(pr3);
  fields_device_pack_yz<B, UNPACK>(mflds, mb, me);
  prof_stop(pr3);
}

template<int B, bool pack>
static void fields_device_pack3_yz(struct psc_mfields *mflds, int mb, int me);

template<int B>
static void
__fields_cuda_from_device3_yz(struct psc_mfields *mflds, int mb, int me)
{
  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_device_pack 3i", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy 3i", 1., 0, 0);
    pr3 = prof_register("field_host_unpack 3i", 1., 0, 0);
  }
    
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(mflds_cuda->ib[1] == -BND);
  assert(mflds_cuda->im[1] >= 2 * B);
  assert(mflds_cuda->im[2] >= 2 * B);

  int nr_map;
  if (B == 4) {
    nr_map = mflds_cuda->nr_map_in;
  } else {
    assert(0);
  }

  prof_start(pr1);
  fields_device_pack3_yz<B, PACK>(mflds, mb, me);
  prof_stop(pr1);
  
  prof_start(pr2);
  assert(B == 4);
  check(cudaMemcpy(mflds_cuda->h_bnd_buf, mflds_cuda->d_bnd_buf,
		   (me - mb) * nr_map * sizeof(*mflds_cuda->h_bnd_buf),
		   cudaMemcpyDeviceToHost));
  prof_stop(pr2);

  prof_start(pr3);
  fields_host_pack3_yz<B, UNPACK>(mflds, mb, me);
  prof_stop(pr3);
}

template<int B>
static void
__fields_cuda_to_device3_yz(struct psc_mfields *mflds, int mb, int me)
{
  static int pr1, pr2, pr3;
  if (!pr1) {
    pr1 = prof_register("field_host_pack 3o", 1., 0, 0);
    pr2 = prof_register("cuda_memcpy 3o", 1., 0, 0);
    pr3 = prof_register("field_device_unpack 3o", 1., 0, 0);
  }
  
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  assert(me - mb <= MAX_BND_COMPONENTS);
  assert(mflds_cuda->ib[1] == -BND);
  assert(mflds_cuda->im[1] >= 2 * B);
  assert(mflds_cuda->im[2] >= 2 * B);

  int nr_map;
  if (B == 2) {
    nr_map = mflds_cuda->nr_map_out;
  } else {
    assert(0);
  }

  prof_start(pr1);
  fields_host_pack3_yz<B, PACK>(mflds, mb, me);
  prof_stop(pr1);

  prof_start(pr2);
  check(cudaMemcpy(mflds_cuda->d_bnd_buf, mflds_cuda->h_bnd_buf,
		   (me - mb) * nr_map * sizeof(*mflds_cuda->d_bnd_buf),
		   cudaMemcpyHostToDevice));
  prof_stop(pr2);

  prof_start(pr3);
  fields_device_pack3_yz<B, UNPACK>(mflds, mb, me);
  prof_stop(pr3);
}

// ======================================================================

EXTERN_C void
__fields_cuda_from_device_inside(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  if (mflds_cuda->im[0] == 2 * -mflds_cuda->ib[0] + 1) {
    __fields_cuda_from_device_yz<2*BND>(mflds, mb, me);
  } else {
    for (int p = 0; p < mflds->nr_patches; p++) {
      struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
      struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
      unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
      check(cudaMemcpy(flds_cuda->bnd.arr,
		       flds_cuda->d_flds + mb * size,
		       (me - mb) * size * sizeof(float),
		       cudaMemcpyDeviceToHost));
    }
  }
}

EXTERN_C void
__fields_cuda_from_device_inside_only(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  if (mflds_cuda->im[0] == 2 * -mflds_cuda->ib[0] + 1) {
    __fields_cuda_from_device3_yz<2*BND>(mflds, mb, me);
  } else {
    for (int p = 0; p < mflds->nr_patches; p++) {
      struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
      struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
      unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
      check(cudaMemcpy(flds_cuda->bnd.arr,
		       flds_cuda->d_flds + mb * size,
		       (me - mb) * size * sizeof(float),
		       cudaMemcpyDeviceToHost));
    }
  }
}

EXTERN_C void
__fields_cuda_to_device_outside(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  if (mflds_cuda->im[0] == 2 * -mflds_cuda->ib[0] + 1) {
    __fields_cuda_to_device3_yz<BND>(mflds, mb, me);
  } else {
    for (int p = 0; p < mflds->nr_patches; p++) {
      struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
      struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
      unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
      check(cudaMemcpy(flds_cuda->d_flds + mb * size,
		       flds_cuda->bnd.arr,
		       (me - mb) * size * sizeof(float),
		       cudaMemcpyHostToDevice));
    }
  }
}

EXTERN_C void
__fields_cuda_to_device_inside(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  if (mflds_cuda->im[0] == 2 * -mflds_cuda->ib[0] + 1) {
    __fields_cuda_to_device_yz<2*BND>(mflds, mb, me);
  } else {
    for (int p = 0; p < mflds->nr_patches; p++) {
      struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
      struct psc_fields_cuda *flds_cuda = psc_fields_cuda(flds);
      unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
      check(cudaMemcpy(flds_cuda->d_flds + mb * size,
		       flds_cuda->bnd.arr,
		       (me - mb) * size * sizeof(float),
		       cudaMemcpyHostToDevice));
    }
  }
}

void
cuda_fill_ghosts_local_gold(float *h_flds, int *nei_patch_by_dir1, int mb, int me,
			    int *im, int nr_fields, int nr_patches, int nr_ghosts, int threadIdx)
{
  int iy, iz;
  int diry, dirz;

  int r_p = threadIdx / nr_ghosts;
  if (r_p >= nr_patches)
    return;
  int tid = threadIdx - nr_ghosts * r_p;

  if (tid < 4 * im[1]) {
    iy = tid % im[1];
    iz = tid / im[1];
    if (iy < 2) {
      diry = -1;
    } else if (iy < im[1] - 2) {
	diry = 0;
    } else {
      diry = 1;
    }
    if (iz < 2) {
      dirz = -1;
    } else {
      dirz = 1;
      iz += im[2] - 2*2;
    }
  } else {
    int tid2 = tid - 4 * im[1];
    iy = tid2 % 4;
    iz = tid2 / 4 + 2;
    dirz = 0;
    if (iy < 2) {
      diry = -1;
    } else {
      diry = 1;
      iy += im[1] - 2*2;
    }
  }
  int s_p = nei_patch_by_dir1[r_p*9 + 3*dirz + diry + 4];
  if (s_p >= 0) {
    float *r_f = &h_flds[((r_p * nr_fields + mb) * im[2]) * im[1]];
    float *s_f = &h_flds[((s_p * nr_fields + mb)
			  * im[2] - dirz * (im[2] - 2*2)) 
			 * im[1] - diry * (im[1] - 2*2)];
    for (int m = 0; m < me - mb; m++) {
      int i = (m * im[2] + iz) * im[1] + iy;
      r_f[i] = s_f[i];
    }
  }
}

static void fields_create_map_out_yz(struct psc_mfields *mflds, int B, int *p_nr_map, int **p_h_map);
static void fields_create_map_in_yz(struct psc_mfields *mflds, int B, int *p_nr_map, int **p_h_map);

EXTERN_C void
__fields_cuda_fill_ghosts_setup(struct psc_mfields *mflds, struct mrc_ddc *ddc)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  int *im = mflds_cuda->im;
  assert(im[0] == 1);
  int nr_patches = mflds->nr_patches;

  if (!mflds_cuda->h_nei_patch) {
    struct mrc_ddc_multi *multi = mrc_ddc_multi(ddc);
    struct mrc_ddc_pattern2 *patt2 = &multi->fill_ghosts2;
    struct mrc_ddc_rank_info *ri = patt2->ri;

    mflds_cuda->h_nei_patch = new int[9 * nr_patches];

    for (int p = 0; p < nr_patches; p++) {
      for (int dir1 = 0; dir1 < 9; dir1++) {
	mflds_cuda->h_nei_patch[p * 9 + dir1] = -1;
      }
    }
    for (int i = 0; i < ri[multi->mpi_rank].n_recv_entries; i++) {
      struct mrc_ddc_sendrecv_entry *re = &ri[multi->mpi_rank].recv_entry[i];
      mflds_cuda->h_nei_patch[re->patch * 9 + re->dir1 / 3] = re->nei_patch;
    }

    check(cudaMalloc((void **) &mflds_cuda->d_nei_patch,
		     9 * nr_patches * sizeof(*mflds_cuda->d_nei_patch)));
    check(cudaMemcpy(mflds_cuda->d_nei_patch, mflds_cuda->h_nei_patch, 
		     9 * nr_patches * sizeof(*mflds_cuda->d_nei_patch),
		     cudaMemcpyHostToDevice));

    fields_create_map_out_yz(mflds, 2, &mflds_cuda->nr_map_out, &mflds_cuda->h_map_out);

    check(cudaMalloc((void **) &mflds_cuda->d_map_out,
		     mflds_cuda->nr_map_out * sizeof(*mflds_cuda->d_map_out)));
    check(cudaMemcpy(mflds_cuda->d_map_out, mflds_cuda->h_map_out, 
		     mflds_cuda->nr_map_out * sizeof(*mflds_cuda->d_map_out),
		     cudaMemcpyHostToDevice));

    fields_create_map_in_yz(mflds, 2, &mflds_cuda->nr_map_in, &mflds_cuda->h_map_in);
    mprintf("map_out %d\n", mflds_cuda->nr_map_out);
    mprintf("map_in %d\n", mflds_cuda->nr_map_in);

    check(cudaMalloc((void **) &mflds_cuda->d_map_in,
		     mflds_cuda->nr_map_in * sizeof(*mflds_cuda->d_map_in)));
    check(cudaMemcpy(mflds_cuda->d_map_in, mflds_cuda->h_map_in, 
		     mflds_cuda->nr_map_in * sizeof(*mflds_cuda->d_map_in),
		     cudaMemcpyHostToDevice));
  }
}

template<int WHAT>
static void
g_fields_device_pack3_yz(int tid, real *d_buf, real *d_flds, int *d_map, int *d_nei_patch_by_dir1,
			 int gmy, int gmz, int nr_patches, int nr_components, int nr_map)
{
  if (tid >= nr_map)
    return;

  // copy only ghost areas that interface with remote patches
  int i = d_map[tid];
  for (int m = 0; m < nr_components; m++) {
    // FIXME, should use F3_DEV_YZ
    if (WHAT == PACK) {
      d_buf[m * nr_map + tid] = d_flds[i + m * gmz * gmy];
    } else if (WHAT == UNPACK) {
      d_flds[i + m * gmz * gmy] = d_buf[m * nr_map + tid]; 
    }
  }
}

template<int WHAT>
__global__ static void
k_fields_device_pack3_yz(real *d_buf, real *d_flds, int *d_map, int *d_nei_patch_by_dir1,
			 int gmy, int gmz, int nr_patches, int nr_components, int nr_map)
{
  int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  if (tid >= nr_map)
    return;

  // copy only ghost areas that interface with remote patches
  int i = d_map[tid];
  for (int m = 0; m < nr_components; m++) {
    // FIXME, should use F3_DEV_YZ
    if (WHAT == PACK) {
      d_buf[m * nr_map + tid] = d_flds[i + m * gmz * gmy];
    } else if (WHAT == UNPACK) {
      d_flds[i + m * gmz * gmy] = d_buf[m * nr_map + tid]; 
    }
  }
}

#undef check

#undef _GLIBCXX_USE_INT128

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

static void
fields_create_map_out_yz(struct psc_mfields *mflds, int B, int *p_nr_map, int **p_h_map)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  bool remote_only = true;
  int *im = mflds_cuda->im;
  int nr_patches = mflds->nr_patches;
  int nr_fields = mflds->nr_fields;

  int nr_map = 0;
  for (int p = 0; p < nr_patches; p++) {
    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	if (dir[1] == 0 && dir[2] == 0) {
	  continue;
	}
	int s_p = mflds_cuda->h_nei_patch[p*9 + 3*dir[2] + dir[1] + 4];
	if (!remote_only || s_p < 0) {
	  nr_map += ((dir[1] == 0 ? im[1] - 2*B : B) *
		     (dir[2] == 0 ? im[2] - 2*B : B));
	}
      }
    }
  }
  *p_nr_map = nr_map;
  *p_h_map = new int[nr_map];

  int tid = 0;
  for (int p = 0; p < nr_patches; p++) {
    for (int jjz = 0; jjz < 2*B; jjz++) {
      int jz = jjz;
      int dirz;
      if (jz < B) {
	dirz = -1;
      } else {
	dirz = 1;
	jz += im[2] - 2*B;
      }
      for (int jy = 0; jy < im[1]; jy++) {
	int diry;
	if (jy < B) {
	  diry = -1;
	} else if (jy < im[1] - B) {
	  diry = 0;
	} else {
	  diry = 1;
	}
	int s_p = mflds_cuda->h_nei_patch[p*9 + 3*dirz + diry + 4];
	if (!remote_only || s_p < 0) {
	  (*p_h_map)[tid++] = ((p * nr_fields + 0) * im[2] + jz) * im[1] + jy;
	}
      }
    }
    for (int jz = B; jz < im[2] - B; jz++) {
      int dirz = 0;
      for (int jjy = 0; jjy < 2*B; jjy++) {
	int jy = jjy;
	int diry;
	if (jy < B) {
	  diry = -1;
	} else {
	  diry = 1;
	  jy += im[1] - 2*B;
	}
	int s_p = mflds_cuda->h_nei_patch[p*9 + 3*dirz + diry + 4];
	if (!remote_only || s_p < 0) {
	  (*p_h_map)[tid++] = ((p * nr_fields + 0) * im[2] + jz) * im[1] + jy;
	}
      }
    }
  }
  assert(tid == nr_map);
}

static void
fields_create_map_in_yz(struct psc_mfields *mflds, int B, int *p_nr_map, int **p_h_map)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  bool remote_only = true;
  int *im = mflds_cuda->im;
  int ldims[3] = { 1, im[1] - 2*2, im[2] - 2*2 };
  int nr_patches = mflds->nr_patches;
  int nr_fields = mflds->nr_fields;
  bool *has_nei = new bool[9 * nr_patches];
  
  // FIXME, we don't need the ghosts here...

  for (int p = 0; p < nr_patches; p++) {
    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	int dir1 = p*9 + 3*dir[2] + dir[1] + 4;
	has_nei[dir1] = false;
	int p_mm = mflds_cuda->h_nei_patch[p*9 + -1 + 3*-1 + 4];
	int p_0m = mflds_cuda->h_nei_patch[p*9 +  0 + 3*-1 + 4];
	int p_pm = mflds_cuda->h_nei_patch[p*9 + +1 + 3*-1 + 4];
	int p_m0 = mflds_cuda->h_nei_patch[p*9 + -1 + 3* 0 + 4];
	int p_p0 = mflds_cuda->h_nei_patch[p*9 + +1 + 3* 0 + 4];
	int p_mp = mflds_cuda->h_nei_patch[p*9 + +1 + 3*+1 + 4];
	int p_0p = mflds_cuda->h_nei_patch[p*9 +  0 + 3*+1 + 4];
	int p_pp = mflds_cuda->h_nei_patch[p*9 + -1 + 3*+1 + 4];
	if (dir[1] == 0 && dir[2] == 0) {
	} else if (dir[1] == -1 && dir[2] == -1) {
	  if (p_mm < 0 || p_m0 < 0 || p_0m < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] ==  0 && dir[2] == -1) {
	  if (p_0m < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == +1 && dir[2] == -1) {
	  if (p_pm < 0 || p_0m < 0 || p_p0 < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == -1 && dir[2] ==  0) {
	  if (p_m0 < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == +1 && dir[2] ==  0) {
	  if (p_p0 < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == -1 && dir[2] == +1) {
	  if (p_mp < 0 || p_m0 < 0 || p_0p < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] ==  0 && dir[2] == +1) {
	  if (p_0p < 0) {
	    has_nei[dir1] = true;
	  }
	} else if (dir[1] == +1 && dir[2] == +1) {
	  if (p_pp < 0 || p_0p < 0 || p_p0 < 0) {
	    has_nei[dir1] = true;
	  }
	}
      }
    }
  }

  int nr_map = 0;
  for (int p = 0; p < nr_patches; p++) {
    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	if (dir[1] == 0 && dir[2] == 0) {
	  continue;
	}
	int dir1 = p*9 + 3*dir[2] + dir[1] + 4;
	if (!remote_only || has_nei[dir1]) {
	  nr_map += ((dir[1] == 0 ? ldims[1] - 2*B : B) *
		     (dir[2] == 0 ? ldims[2] - 2*B : B));
	}
      }
    }
  }
  *p_nr_map = nr_map;
  *p_h_map = new int[nr_map];

  int tid = 0;
  for (int p = 0; p < nr_patches; p++) {
    for (int jjz = 0; jjz < 2*B; jjz++) {
      int jz = jjz;
      int dirz;
      if (jz < B) {
	dirz = -1;
      } else {
	dirz = 1;
	jz += ldims[2] - 2*B;
      }
      for (int jy = 0; jy < ldims[1]; jy++) {
	int diry;
	if (jy < B) {
	  diry = -1;
	} else if (jy < ldims[1] - B) {
	  diry = 0;
	} else {
	  diry = 1;
	}
	int dir1 = p*9 + 3*dirz + diry + 4;
	if (!remote_only || has_nei[dir1]) {
	  (*p_h_map)[tid++] = ((p * nr_fields + 0) * im[2] + (jz + B)) * im[1] + (jy + B);
	}
      }
    }
    for (int jz = B; jz < ldims[2] - B; jz++) {
      int dirz = 0;
      for (int jjy = 0; jjy < 2*B; jjy++) {
	int jy = jjy;
	int diry;
	if (jy < B) {
	  diry = -1;
	} else {
	  diry = 1;
	  jy += ldims[1] - 2*B;
	}
	int dir1 = p*9 + 3*dirz + diry + 4;
	if (!remote_only || has_nei[dir1]) {
	  (*p_h_map)[tid++] = ((p * nr_fields + 0) * im[2] + (jz + B)) * im[1] + (jy + B);
	}
      }
    }
  }
  assert(tid == nr_map);
  delete[] has_nei;
}

template<int B, bool pack>
static void
fields_device_pack3_yz(struct psc_mfields *mflds, int mb, int me)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);

  int *im = mflds_cuda->im;
  assert(im[0] == 1);
  int nr_patches = mflds->nr_patches;
  int nr_map;
  int *d_map;
  if (B == 2) {
    nr_map = mflds_cuda->nr_map_out;
    d_map = mflds_cuda->d_map_out;
  } else if (B == 4) {
    nr_map = mflds_cuda->nr_map_in;
    d_map = mflds_cuda->d_map_in;
  } else {
    assert(0);
  }

#if 1
  dim3 dimGrid((nr_map + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);
    
  float *d_flds = mflds_cuda->d_flds + mb * im[1] * im[2];
  if (me - mb == 3) {
    k_fields_device_pack3_yz<pack> <<<dimGrid, dimBlock>>>
      (mflds_cuda->d_bnd_buf, d_flds, d_map, mflds_cuda->d_nei_patch,
       im[1], im[2], nr_patches, 3, nr_map);
  } else if (me - mb == 1) {
    k_fields_device_pack3_yz<pack> <<<dimGrid, dimBlock>>>
      (mflds_cuda->d_bnd_buf, d_flds, d_map, mflds_cuda->d_nei_patch,
       im[1], im[2], nr_patches, 1, nr_map);
  } else {
    assert(0);
  }
  cuda_sync_if_enabled();
#else
  thrust::device_ptr<float> d_bnd_buf(mflds_cuda->d_bnd_buf);
  thrust::device_ptr<float> d_flds(mflds_cuda->d_flds);
  thrust::host_vector<float> h_bnd_buf(d_bnd_buf, d_bnd_buf + nr_patches * nr_ghosts * NR_COMPONENTS);
  thrust::host_vector<float> h_flds(d_flds, d_flds + nr_patches * nr_fields * im[1] * im[2]);

  for (int tid = 0; tid < nr_map; tid++) {
    g_fields_device_pack3_yz<pack>
      (tid, &h_bnd_buf[0], &h_flds[mb * im[1] * im[2]], &h_map[0], &h_nei_patch_by_dir1[0],
       im[1], im[2], nr_patches, NR_COMPONENTS, nr_map);
  }

  thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
#endif
}

