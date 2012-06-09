
#include <psc_cuda.h>
#include "cuda_sort2.h"
#include "particles_cuda.h"
#include "psc_bnd_cuda.h"

#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>

#define PFX(x) xchg_##x
#include "constants.c"

// FIXME const mem for dims?
// FIXME probably should do our own loop rather than use blockIdx

__global__ static void
exchange_particles(int n_part, particles_cuda_dev_t h_dev,
		   int ldimsx, int ldimsy, int ldimsz)
{
  int ldims[3] = { ldimsx, ldimsy, ldimsz };
  int xm[3];

  for (int d = 0; d < 3; d++) {
    xm[d] = ldims[d] / d_consts.dxi[d];
  }

  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (i < n_part) {
    particle_cuda_real_t xi[3] = {
      h_dev.xi4[i].x * d_consts.dxi[0],
      h_dev.xi4[i].y * d_consts.dxi[1],
      h_dev.xi4[i].z * d_consts.dxi[2] };
    int pos[3];
    for (int d = 0; d < 3; d++) {
      pos[d] = cuda_fint(xi[d]);
    }
    if (pos[1] < 0) {
      h_dev.xi4[i].y += xm[1];
      if (h_dev.xi4[i].y >= xm[1])
	h_dev.xi4[i].y = 0.f;
    }
    if (pos[2] < 0) {
      h_dev.xi4[i].z += xm[2];
      if (h_dev.xi4[i].z >= xm[2])
	h_dev.xi4[i].z = 0.f;
    }
    if (pos[1] >= ldims[1]) {
      h_dev.xi4[i].y -= xm[1];
    }
    if (pos[2] >= ldims[2]) {
      h_dev.xi4[i].z -= xm[2];
    }
  }
}

EXTERN_C void
cuda_exchange_particles(int p, struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  struct psc_patch *patch = &ppsc->patch[p];

  xchg_set_constants(prts, NULL);

  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (prts->n_part + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     exchange_particles, (prts->n_part, *cuda->h_dev,
				  patch->ldims[0], patch->ldims[1], patch->ldims[2]));
}

// ----------------------------------------------------------------------
// cuda_mprts_find_block_indices_2_total
//
// like cuda_find_block_indices, but handles out-of-bound
// particles

__global__ static void
mprts_find_block_indices_2_total(struct cuda_params prm, particles_cuda_dev_t *d_cp_prts,
				 unsigned int *d_bidx, int nr_patches)
{
  int tid = threadIdx.x;

  int block_pos[3];
  block_pos[1] = blockIdx.x;
  block_pos[2] = blockIdx.y % prm.b_mx[2];
  int bid = block_pos_to_block_idx(block_pos, prm.b_mx);
  int p = blockIdx.y / prm.b_mx[2];

  // FIXME/OPT, could be done better like reorder_send_buf
  int block_begin = d_cp_prts[p].d_off[bid];
  int block_end   = d_cp_prts[p].d_off[bid+1];

  int nr_blocks = prm.b_mx[1] * prm.b_mx[2];

  for (int n = block_begin + tid; n < block_end; n += THREADS_PER_BLOCK) {
    float4 xi4 = d_cp_prts[0].xi4[n];
    unsigned int block_pos_y = cuda_fint(xi4.y * prm.b_dxi[1]);
    unsigned int block_pos_z = cuda_fint(xi4.z * prm.b_dxi[2]);

    int block_idx;
    if (block_pos_y >= prm.b_mx[1] || block_pos_z >= prm.b_mx[2]) {
      block_idx = nr_blocks * nr_patches;
    } else {
      block_idx = block_pos_z * prm.b_mx[1] + block_pos_y + p * nr_blocks;
    }
    d_bidx[n] = block_idx;
  }
}

EXTERN_C void
cuda_mprts_find_block_indices_2_total(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  if (mprts->nr_patches == 0) {
    return;
  }

  struct cuda_params prm;
  set_params(&prm, ppsc, psc_mparticles_get_patch(mprts, 0), NULL);
    
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { prm.b_mx[1], prm.b_mx[2] * mprts->nr_patches };
  
  RUN_KERNEL(dimGrid, dimBlock,
	     mprts_find_block_indices_2_total, (prm, psc_mparticles_cuda(mprts)->d_dev,
						mprts_cuda->d_bidx, mprts->nr_patches));
  free_params(&prm);
}

// ----------------------------------------------------------------------
// cuda_mprts_find_block_indices_ids_total

__global__ static void
mprts_find_block_indices_ids_total(struct cuda_params prm, particles_cuda_dev_t *d_cp_prts,
				   unsigned int *d_bidx, unsigned int *d_ids, int nr_patches)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  int nr_blocks = prm.b_mx[1] * prm.b_mx[2];

  unsigned int off = 0;
  for (int p = 0; p < nr_patches; p++) {
    if (n < d_cp_prts[p].n_part) {
      float4 xi4 = d_cp_prts[0].xi4[n + off];
      unsigned int block_pos_y = cuda_fint(xi4.y * prm.b_dxi[1]);
      unsigned int block_pos_z = cuda_fint(xi4.z * prm.b_dxi[2]);
      
      int block_idx;
      if (block_pos_y >= prm.b_mx[1] || block_pos_z >= prm.b_mx[2]) {
	block_idx = -1; // not supposed to happen here!
      } else {
	block_idx = block_pos_z * prm.b_mx[1] + block_pos_y + p * nr_blocks;
      }
      d_bidx[n + off] = block_idx;
      d_ids[n + off] = n + off;
    }
    off += d_cp_prts[p].n_part;
  }
}

EXTERN_C void
cuda_mprts_find_block_indices_ids_total(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  if (mprts->nr_patches == 0) {
    return;
  }

  int max_n_part = 0;
  mprts_cuda->nr_prts_send = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    mprts_cuda->nr_prts_send += cuda->bnd_n_send;
    if (prts->n_part > max_n_part) {
      max_n_part = prts->n_part;
    }
  }

  struct cuda_params prm;
  set_params(&prm, ppsc, psc_mparticles_get_patch(mprts, 0), NULL);
    
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (max_n_part + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };

  RUN_KERNEL(dimGrid, dimBlock,
	     mprts_find_block_indices_ids_total, (prm, psc_mparticles_cuda(mprts)->d_dev,
						  mprts_cuda->d_bidx, mprts_cuda->d_ids,
						  mprts->nr_patches));
  free_params(&prm);
}

// ======================================================================
// cuda_mprts_find_block_indices_3

EXTERN_C void
cuda_mprts_find_block_indices_3(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);

    off += prts->n_part;
  }
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    // for consistency, use same block indices that we counted earlier
    check(cudaMemcpy(mprts_cuda->d_bidx + off, cuda->bnd_idx,
		     cuda->bnd_n_recv * sizeof(*mprts_cuda->d_bidx),
		     cudaMemcpyHostToDevice));
    // abuse of alt_bidx!!! FIXME
    check(cudaMemcpy(mprts_cuda->d_alt_bidx + off, cuda->bnd_off,
		     cuda->bnd_n_recv * sizeof(*mprts_cuda->d_alt_bidx),
		     cudaMemcpyHostToDevice));
    off += cuda->bnd_n_recv;
  }
  assert(off == mprts_cuda->nr_prts);
}

// ======================================================================
// mprts_reorder_send_buf_total

__global__ static void
mprts_reorder_send_buf_total(int nr_prts, int nr_oob, unsigned int *d_bidx, unsigned int *d_sums,
			     float4 *d_xi4, float4 *d_pxi4,
			     float4 *d_xchg_xi4, float4 *d_xchg_pxi4)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (i >= nr_prts)
    return;

  if (d_bidx[i] == nr_oob) {
    int j = d_sums[i];
    d_xchg_xi4[j]  = d_xi4[i];
    d_xchg_pxi4[j] = d_pxi4[i];
  }
}

EXTERN_C void
cuda_mprts_reorder_send_buf_total(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  if (mprts->nr_patches == 0)
    return;

  struct psc_particles_cuda *cuda = psc_particles_cuda(psc_mparticles_get_patch(mprts, 0));
  
  mprts_cuda->nr_prts_send = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    mprts_cuda->nr_prts_send += cuda->bnd_n_send;
  }

  float4 *xchg_xi4 = mprts_cuda->d_xi4 + mprts_cuda->nr_prts;
  float4 *xchg_pxi4 = mprts_cuda->d_pxi4 + mprts_cuda->nr_prts;
  assert(mprts_cuda->nr_prts + mprts_cuda->nr_prts_send < mprts_cuda->nr_alloced);
  int nr_oob = cuda->nr_blocks * mprts->nr_patches;
  
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (mprts_cuda->nr_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };
  
  RUN_KERNEL(dimGrid, dimBlock,
	     mprts_reorder_send_buf_total, (mprts_cuda->nr_prts, nr_oob,
					    mprts_cuda->d_bidx, mprts_cuda->d_sums,
					    mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
					    xchg_xi4, xchg_pxi4));
}

// ======================================================================
// psc_mparticles_cuda_swap_alt

static void
psc_mparticles_cuda_swap_alt(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    float4 *alt_xi4 = cuda->h_dev->alt_xi4;
    float4 *alt_pxi4 = cuda->h_dev->alt_pxi4;
    cuda->h_dev->alt_xi4 = cuda->h_dev->xi4;
    cuda->h_dev->alt_pxi4 = cuda->h_dev->pxi4;
    cuda->h_dev->xi4 = alt_xi4;
    cuda->h_dev->pxi4 = alt_pxi4;
  }
  float4 *tmp_xi4 = mprts_cuda->d_alt_xi4;
  float4 *tmp_pxi4 = mprts_cuda->d_alt_pxi4;
  mprts_cuda->d_alt_xi4 = mprts_cuda->d_xi4;
  mprts_cuda->d_alt_pxi4 = mprts_cuda->d_pxi4;
  mprts_cuda->d_xi4 = tmp_xi4;
  mprts_cuda->d_pxi4 = tmp_pxi4;
}

// ======================================================================
// reorder_and_offsets

__global__ static void
mprts_reorder_and_offsets(int nr_prts, float4 *xi4, float4 *pxi4, float4 *alt_xi4, float4 *alt_pxi4,
			  unsigned int *d_bidx, unsigned int *d_ids, unsigned int *d_off, int last_block)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (i > nr_prts)
    return;

  int block, prev_block;
  if (i < nr_prts) {
    alt_xi4[i] = xi4[d_ids[i]];
    alt_pxi4[i] = pxi4[d_ids[i]];
    
    block = d_bidx[i];
  } else { // needed if there is no particle in the last block
    block = last_block;
  }

  // OPT: d_bidx[i-1] could use shmem
  // create offsets per block into particle array
  prev_block = -1;
  if (i > 0) {
    prev_block = d_bidx[i-1];
  }
  for (int b = prev_block + 1; b <= block; b++) {
    d_off[b] = i;
  }
}

void
cuda_mprts_reorder_and_offsets(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  if (mprts->nr_patches == 0) {
    return;
  }
  int nr_blocks = psc_particles_cuda(psc_mparticles_get_patch(mprts, 0))->nr_blocks;

  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (mprts_cuda->nr_prts + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     mprts_reorder_and_offsets, (mprts_cuda->nr_prts, mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
					 mprts_cuda->d_alt_xi4, mprts_cuda->d_alt_pxi4,
					 mprts_cuda->d_bidx, mprts_cuda->d_ids,
					 mprts_cuda->d_off, mprts->nr_patches * nr_blocks));

  psc_mparticles_cuda_swap_alt(mprts);
  psc_mparticles_cuda_copy_to_dev(mprts);
}

// ======================================================================
// cuda_mprts_copy_from_dev

void
cuda_mprts_copy_from_dev(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  if (mprts->nr_patches == 0) {
    return;
  }

  mprts_cuda->h_bnd_xi4 = new float4[mprts_cuda->nr_prts_send];
  mprts_cuda->h_bnd_pxi4 = new float4[mprts_cuda->nr_prts_send];

  check(cudaMemcpy(mprts_cuda->h_bnd_xi4, mprts_cuda->d_xi4 + mprts_cuda->nr_prts,
		   mprts_cuda->nr_prts_send * sizeof(float4), cudaMemcpyDeviceToHost));
  check(cudaMemcpy(mprts_cuda->h_bnd_pxi4, mprts_cuda->d_pxi4 + mprts_cuda->nr_prts,
		   mprts_cuda->nr_prts_send * sizeof(float4), cudaMemcpyDeviceToHost));
}

//======================================================================
// cuda_mprts_convert_from_cuda

void
cuda_mprts_convert_from_cuda(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  if (mprts->nr_patches == 0) {
    return;
  }

  float4 *bnd_xi4 = mprts_cuda->h_bnd_xi4;
  float4 *bnd_pxi4 = mprts_cuda->h_bnd_pxi4;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    cuda->bnd_prts = new particle_single_t[cuda->bnd_n_send];
    for (int n = 0; n < cuda->bnd_n_send; n++) {
      particle_single_t *prt = &cuda->bnd_prts[n];
      prt->xi      = bnd_xi4[n].x;
      prt->yi      = bnd_xi4[n].y;
      prt->zi      = bnd_xi4[n].z;
      prt->kind    = cuda_float_as_int(bnd_xi4[n].w);
      prt->pxi     = bnd_pxi4[n].x;
      prt->pyi     = bnd_pxi4[n].y;
      prt->pzi     = bnd_pxi4[n].z;
      prt->qni_wni = bnd_pxi4[n].w;
    }
    bnd_xi4 += cuda->bnd_n_send;
    bnd_pxi4 += cuda->bnd_n_send;
  }
  delete[] mprts_cuda->h_bnd_xi4;
  delete[] mprts_cuda->h_bnd_pxi4;
}

// ======================================================================
// cuda_mprts_copy_to_dev

void
cuda_mprts_copy_to_dev(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  float4 *d_xi4 = mprts_cuda->d_xi4;
  float4 *d_pxi4 = mprts_cuda->d_pxi4;

  unsigned int off = mprts_cuda->nr_prts;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    assert(off + cuda->bnd_n_recv <= mprts_cuda->nr_alloced);

    check(cudaMemcpy(d_xi4 + off, cuda->bnd_xi4,
		     cuda->bnd_n_recv * sizeof(*cuda->bnd_xi4),
		     cudaMemcpyHostToDevice));
    check(cudaMemcpy(d_pxi4 + off, cuda->bnd_pxi4,
		     cuda->bnd_n_recv * sizeof(*cuda->bnd_pxi4),
		     cudaMemcpyHostToDevice));

    off += cuda->bnd_n_recv;
  }
  mprts_cuda->nr_prts = off;
}

// ======================================================================
// cuda_mprts_sort

void
cuda_mprts_sort(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  unsigned int *h_bidx = new unsigned int[mprts_cuda->nr_alloced];
  check(cudaMemcpy(h_bidx, mprts_cuda->d_bidx, mprts_cuda->nr_prts * sizeof(unsigned int),
		   cudaMemcpyDeviceToHost));

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    off += prts->n_part;
  }
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    for (int n = 0; n < cuda->bnd_n_recv; n++) {
      assert(h_bidx[off + n] < cuda->nr_blocks);
      h_bidx[off + n] += p * cuda->nr_blocks;
    }
    off += cuda->bnd_n_recv;
  }
  assert(off = mprts_cuda->nr_prts);

  check(cudaMemcpy(mprts_cuda->d_bidx, h_bidx, mprts_cuda->nr_prts * sizeof(unsigned int),
		   cudaMemcpyHostToDevice));

  cuda_mprts_sort_pairs_device(mprts);
#if 0
    // OPT: when calculating bidx, do preprocess then
    void *sp = sort_pairs_3_create(cuda->b_mx);
    sort_pairs_3_device(sp, mprts_cuda->d_bidx + off, mprts_cuda->d_alt_bidx + off,
			mprts_cuda->d_ids + off, prts->n_part, cuda->h_dev->offsets,
			cuda->bnd_n_part_save, cuda->bnd_cnt);
    sort_pairs_3_destroy(sp);
#endif

  unsigned int *h_ids = new unsigned int[mprts_cuda->nr_alloced];
  check(cudaMemcpy(h_ids, mprts_cuda->d_ids, mprts_cuda->nr_alloced * sizeof(unsigned int),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(h_bidx, mprts_cuda->d_bidx, mprts_cuda->nr_alloced * sizeof(unsigned int),
		   cudaMemcpyDeviceToHost));

  off = 0;
  int last = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    prts->n_part += cuda->bnd_n_recv - cuda->bnd_n_send;

    for (int i = 0; i < prts->n_part; i++) {
      if (h_bidx[off + i] < last) {
	mprintf("i %d last %d bidx %d\n", i, last, h_bidx[off + i]);
      }
      assert(h_bidx[off + i] >= last);
      assert(h_bidx[off + i] >= p * cuda->nr_blocks);
      if (!(h_bidx[off + i] < (p + 1) * cuda->nr_blocks)) {
	mprintf("i %d h_bidx %d p %d\n", i, h_bidx[off + i], p);
      }
      assert(h_bidx[off + i] < (p + 1) * cuda->nr_blocks);
      assert(h_ids[off + i] < mprts_cuda->nr_prts);
      last = h_bidx[off + i];
    }

    cuda->h_dev->xi4 = mprts_cuda->d_xi4 + off;
    cuda->h_dev->pxi4 = mprts_cuda->d_pxi4 + off;
    cuda->h_dev->alt_xi4 = mprts_cuda->d_alt_xi4 + off;
    cuda->h_dev->alt_pxi4 = mprts_cuda->d_alt_pxi4 + off;
    cuda->n_alloced = prts->n_part;
    off += prts->n_part;
  }
  mprts_cuda->nr_prts = off;
  delete[] h_bidx;

  psc_mparticles_cuda_copy_to_dev(mprts);
}

// ======================================================================
// cuda_mprts_reorder

__global__ static void
mprts_reorder(int nr_prts, unsigned int *d_ids,
	      float4 *xi4, float4 *pxi4,
	      float4 *alt_xi4, float4 *alt_pxi4)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (i < nr_prts) {
    int j = d_ids[i];
    alt_xi4[i] = xi4[j];
    alt_pxi4[i] = pxi4[j];
  }
}

void
cuda_mprts_reorder(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (mprts_cuda->nr_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     mprts_reorder, (mprts_cuda->nr_prts, mprts_cuda->d_ids,
			     mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
			     mprts_cuda->d_alt_xi4, mprts_cuda->d_alt_pxi4));
  
  psc_mparticles_cuda_swap_alt(mprts);
}

// ======================================================================
// cuda_mprts_free

void
cuda_mprts_free(struct psc_mparticles *mprts)
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    free(cuda->bnd_idx);
    free(cuda->bnd_off);
    free(cuda->bnd_cnt);
    free(cuda->bnd_prts);
    free(cuda->bnd_xi4);
    free(cuda->bnd_pxi4);
  }
}

// ======================================================================
// cuda_mprts_check_ordered_total

void
cuda_mprts_check_ordered_total(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  psc_mparticles_cuda_copy_to_dev(mprts); // update n_part, particle pointers
  cuda_mprts_find_block_indices_2_total(mprts);

  unsigned int last = 0;
  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    unsigned int *bidx = new unsigned int[prts->n_part];
    cuda_copy_bidx_from_dev(prts, bidx, mprts_cuda->d_bidx + off);
    
    for (int n = 0; n < prts->n_part; n++) {
      if (!(bidx[n] >= last && bidx[n] < mprts->nr_patches * cuda->nr_blocks)) {
	mprintf("p = %d, n = %d bidx = %d last = %d\n", p, n, bidx[n], last);
	assert(0);
      }
      last = bidx[n];
    }

    delete[] bidx;

    off += prts->n_part;
  }
}

// ======================================================================
// cuda_mprts_compact

void
cuda_mprts_compact(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  float4 *d_alt_xi4 = mprts_cuda->d_alt_xi4;
  float4 *d_alt_pxi4 = mprts_cuda->d_alt_pxi4;
  float4 *d_xi4 = mprts_cuda->d_xi4;
  float4 *d_pxi4 = mprts_cuda->d_pxi4;

  int nr_prts = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    assert(d_alt_xi4 + prts->n_part <= mprts_cuda->d_alt_xi4 + mprts_cuda->nr_alloced);
    check(cudaMemcpy(d_alt_xi4, cuda->h_dev->xi4,
		     prts->n_part * sizeof(*cuda->h_dev->alt_xi4),
		     cudaMemcpyDeviceToDevice));
    check(cudaMemcpy(d_alt_pxi4, cuda->h_dev->pxi4,
		     prts->n_part * sizeof(*cuda->h_dev->alt_pxi4),
		     cudaMemcpyDeviceToDevice));
    nr_prts += prts->n_part;
    cuda->n_alloced = prts->n_part;
    cuda->h_dev->alt_xi4 = d_alt_xi4;
    cuda->h_dev->alt_pxi4 = d_alt_pxi4;
    cuda->h_dev->xi4 = d_xi4;
    cuda->h_dev->pxi4 = d_pxi4;
    d_alt_xi4 += cuda->n_alloced;
    d_alt_pxi4 += cuda->n_alloced;
    d_xi4 += cuda->n_alloced;
    d_pxi4 += cuda->n_alloced;
  }
  mprts_cuda->nr_prts = nr_prts;
  psc_mparticles_cuda_swap_alt(mprts);
  psc_mparticles_cuda_copy_to_dev(mprts);
}

