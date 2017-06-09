
#include <cstdio>
#include <cassert>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#define cudaCheck(ierr) do {						\
    if (ierr != cudaSuccess)						\
      fprintf(stderr, "IERR = %d (%s)\n", ierr, cudaGetErrorName(ierr)); \
    assert(ierr == cudaSuccess);					\
  } while(0)

struct cuda_domain_info {
  int nr_patches;
  int mx[3]; // number of cells per patch
  int bs[3]; // size of each block / super-cell
  double dx[3]; // size of a single cell
};

void
cuda_domain_info_set_test_1(struct cuda_domain_info *info)
{
  info->nr_patches = 1;
  info->mx[0] = 1; info->mx[1] = 4; info->mx[2] = 2;
  info->bs[0] = 1; info->bs[1] = 1; info->bs[2] = 1;
  info->dx[0] = 1.; info->dx[1] = 10.; info->dx[2] = 10.;
  for (int d = 0; d < 3; d++) {
    assert(info->mx[d] % info->bs[d] == 0);
  }
};

struct cuda_mparticles {
  unsigned int nr_prts;     // total # of particles in all patches
  unsigned int nr_alloced;  // arrays are alloced for this # of particles
  unsigned int nr_patches;
  unsigned int nr_blocks;
  unsigned int *nr_prts_by_patch;
  int mx[3];      // number of cells per direction in each patch
  int b_mx[3];    // number of blocks per direction in each patch
  float dx[3];    // cell size (in actual length units)
  float b_dxi[3]; // inverse of block size (in actual length units)

  // per particle
  float4 *d_xi4, *d_pxi4;
  float4 *d_alt_xi4, *d_alt_pxi4;
  unsigned int *d_bidx;
  unsigned int *d_id;

  // per patch
  unsigned int *d_nr_prts_by_patch;

  // per block
  unsigned int *d_off;
};

void
cuda_mparticles_set_domain_info(struct cuda_mparticles *cuda_mprts,
				const struct cuda_domain_info *info)
{
  cuda_mprts->nr_patches = info->nr_patches;
  for (int d = 0; d < 3; d++) {
    cuda_mprts->mx[d] = info->mx[d];
    cuda_mprts->b_mx[d] = info->mx[d] / info->bs[d];
    cuda_mprts->dx[d] = info->dx[d];
    cuda_mprts->b_dxi[d] = 1.f / (info->bs[d] * info->dx[d]);
  }
  cuda_mprts->nr_blocks = info->nr_patches *
    cuda_mprts->b_mx[0] * cuda_mprts->b_mx[1] * cuda_mprts->b_mx[2];
}

void
cuda_mparticles_alloc(struct cuda_mparticles *cuda_mprts, int nr_prts)
{
  cudaError_t ierr;
  unsigned int nr_alloced = nr_prts * 1.4;
  cuda_mprts->nr_prts = nr_prts;
  cuda_mprts->nr_alloced = nr_alloced;
  cuda_mprts->nr_prts_by_patch = new unsigned int[cuda_mprts->nr_patches];

  assert(cuda_mprts->nr_patches == 1);
  for (int p = 0; p < cuda_mprts->nr_patches; p++) {
    cuda_mprts->nr_prts_by_patch[p] = nr_prts;
  }

  ierr = cudaMalloc((void **) &cuda_mprts->d_xi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_pxi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_alt_xi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_alt_pxi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_bidx, nr_alloced * sizeof(unsigned int)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_id, nr_alloced * sizeof(unsigned int)); cudaCheck(ierr);

  ierr = cudaMalloc((void **) &cuda_mprts->d_nr_prts_by_patch, cuda_mprts->nr_patches * sizeof(unsigned int)); cudaCheck(ierr);

  ierr = cudaMalloc((void **) &cuda_mprts->d_off, (cuda_mprts->nr_blocks + 1) * sizeof(unsigned int)); cudaCheck(ierr);
}

void
cuda_mparticles_free(struct cuda_mparticles *cuda_mprts)
{
  cudaError_t ierr;

  delete[] cuda_mprts->nr_prts_by_patch;
  
  ierr = cudaFree(cuda_mprts->d_xi4); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_pxi4); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_alt_xi4); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_alt_pxi4); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_bidx); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_id); cudaCheck(ierr);

  ierr = cudaFree(cuda_mprts->d_nr_prts_by_patch); cudaCheck(ierr);

  ierr = cudaFree(cuda_mprts->d_off); cudaCheck(ierr);
}

void
cuda_mparticles_copy_to_dev(struct cuda_mparticles *cuda_mprts)
{
  cudaError_t ierr;
  
  ierr = cudaMemcpy(cuda_mprts->d_nr_prts_by_patch, cuda_mprts->nr_prts_by_patch,
		    cuda_mprts->nr_patches * sizeof(unsigned int),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);
}

void
cuda_mparticles_set_test_1(struct cuda_mparticles *cuda_mprts)
{
  int nr_prts = cuda_mprts->mx[0] * cuda_mprts->mx[1] * cuda_mprts->mx[2]
    * cuda_mprts->nr_patches;
  cuda_mparticles_alloc(cuda_mprts, nr_prts);
  
  thrust::device_ptr<float4> d_xi4(cuda_mprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cuda_mprts->d_pxi4);

  int *mx = cuda_mprts->mx;
  float *dx = cuda_mprts->dx;
  
  int n = 0;
  for (int i = 0; i < mx[0]; i++) {
    for (int j = 0; j < mx[1]; j++) {
      for (int k = 0; k < mx[2]; k++) {
	d_xi4[n] = (float4) { dx[0] * (i + .5f),
			      dx[1] * (j + .5f),
			      dx[1] * (k + .5f), 1. };
	d_pxi4[n] = (float4) { i, j, k, 2. };
	n++;
      }
    }
  }
}

void
cuda_mparticles_dump(struct cuda_mparticles *cuda_mprts)
{
  int nr_prts = cuda_mprts->nr_prts;
  
  thrust::device_ptr<float4> d_xi4(cuda_mprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cuda_mprts->d_pxi4);
  thrust::device_ptr<unsigned int> d_bidx(cuda_mprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cuda_mprts->d_id);
  thrust::device_ptr<unsigned int> d_off(cuda_mprts->d_off);

  printf("cuda_mparticles_dump: nr_prts = %d\n", nr_prts); 
  for (int n = 0; n < nr_prts; n++) {
    float4 xi4 = d_xi4[n], pxi4 = d_pxi4[n];
    unsigned int bidx = d_bidx[n], id = d_id[n];
    printf("cuda_mparticles_dump: [%d] %g %g %g // %g // %g %g %g // %g || bidx %d id %d\n",
	   n, xi4.x, xi4.y, xi4.z, xi4.w, pxi4.x, pxi4.y, pxi4.z, pxi4.w,
	   bidx, id);
  }

  for (int b = 0; b <= cuda_mprts->nr_blocks; b++) {
    unsigned int off = d_off[b];
    printf("cuda_mparticles_dump: off[%d] = %d\n", b, off);
  }
}

// ----------------------------------------------------------------------
// cuda_params

struct cuda_params {
  unsigned int b_mx[3];
  float b_dxi[3];
};

void
cuda_params_set(struct cuda_params *prm, const struct cuda_mparticles *cuda_mprts)
{
  for (int d = 0; d < 3; d++) {
    prm->b_mx[d]  = cuda_mprts->b_mx[d];
    prm->b_dxi[d] = cuda_mprts->b_dxi[d];
  }
}

void
cuda_params_free(struct cuda_params *prm)
{
}

// ----------------------------------------------------------------------
// cuda_mparticles_find_block_indices_ids_total

#define THREADS_PER_BLOCK 512

__global__ static void
mprts_find_block_indices_ids_total(struct cuda_params prm, float4 *d_xi4, unsigned int *d_nr_prts,
				   unsigned int *d_bidx, unsigned int *d_id, int nr_patches)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  int nr_blocks = prm.b_mx[1] * prm.b_mx[2];

  unsigned int off = 0;
  for (int p = 0; p < nr_patches; p++) {
    if (n < d_nr_prts[p]) {
      float4 xi4 = d_xi4[n + off];
      unsigned int block_pos_y = __float2int_rd(xi4.y * prm.b_dxi[1]);
      unsigned int block_pos_z = __float2int_rd(xi4.z * prm.b_dxi[2]);
      
      int block_idx;
      if (block_pos_y >= prm.b_mx[1] || block_pos_z >= prm.b_mx[2]) {
	block_idx = -1; // not supposed to happen here!
      } else {
	block_idx = block_pos_z * prm.b_mx[1] + block_pos_y + p * nr_blocks;
      }
      d_bidx[n + off] = block_idx;
      d_id[n + off] = n + off;
    }
    off += d_nr_prts[p];
  }
}

void
cuda_mparticles_find_block_indices_ids_total(struct cuda_mparticles *cuda_mprts)
{
  if (cuda_mprts->nr_patches == 0) {
    return;
  }

  int max_nr_prts = 0;
  int nr_prts = 0;
  //mprts_cuda->nr_prts_send = 0;
  for (int p = 0; p < cuda_mprts->nr_patches; p++) {
    if (cuda_mprts->nr_prts_by_patch[p] > max_nr_prts) {
      max_nr_prts = cuda_mprts->nr_prts_by_patch[p];
    }
    nr_prts += cuda_mprts->nr_prts_by_patch[p];
    // struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    // struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    // mprts_cuda->nr_prts_send += cuda->bnd_n_send;
    // if (prts->n_part > max_n_part) {
    //   max_n_part = prts->n_part;
    // }
    // cuda->h_dev->n_part = prts->n_part;
    // nr_prts += prts->n_part;
  }
  //mprts_cuda->nr_prts = nr_prts;
  assert(cuda_mprts->nr_prts == nr_prts);
  cuda_mparticles_copy_to_dev(cuda_mprts);

  struct cuda_params prm;
  cuda_params_set(&prm, cuda_mprts);
    
  dim3 dimGrid((max_nr_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  mprts_find_block_indices_ids_total<<<dimGrid, dimBlock>>>(prm, cuda_mprts->d_xi4, 
							    cuda_mprts->d_nr_prts_by_patch,
							    cuda_mprts->d_bidx,
							    cuda_mprts->d_id,
							    cuda_mprts->nr_patches);
  cuda_params_free(&prm);
}

// ======================================================================
// cuda_mparticles_swap_alt

void
cuda_mparticles_swap_alt(struct cuda_mparticles *cuda_mprts)
{
  float4 *tmp_xi4 = cuda_mprts->d_alt_xi4;
  float4 *tmp_pxi4 = cuda_mprts->d_alt_pxi4;
  cuda_mprts->d_alt_xi4 = cuda_mprts->d_xi4;
  cuda_mprts->d_alt_pxi4 = cuda_mprts->d_pxi4;
  cuda_mprts->d_xi4 = tmp_xi4;
  cuda_mprts->d_pxi4 = tmp_pxi4;
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
cuda_mparticles_reorder_and_offsets(struct cuda_mparticles *cuda_mprts)
{
  if (cuda_mprts->nr_patches == 0) {
    return;
  }

  dim3 dimGrid((cuda_mprts->nr_prts + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  mprts_reorder_and_offsets<<<dimGrid, dimBlock>>>(cuda_mprts->nr_prts, cuda_mprts->d_xi4, cuda_mprts->d_pxi4,
						   cuda_mprts->d_alt_xi4, cuda_mprts->d_alt_pxi4,
						   cuda_mprts->d_bidx, cuda_mprts->d_id,
						   cuda_mprts->d_off, cuda_mprts->nr_blocks);

  cuda_mparticles_swap_alt(cuda_mprts);
  //  psc_mparticles_cuda_copy_to_dev(mprts);
}



int
main(void)
{
  struct cuda_mparticles _cuda_mprts, *cuda_mprts = &_cuda_mprts;

  struct cuda_domain_info info;
  cuda_domain_info_set_test_1(&info);

  cuda_mparticles_set_domain_info(cuda_mprts, &info);
  cuda_mparticles_set_test_1(cuda_mprts);
  cuda_mparticles_dump(cuda_mprts);

  cuda_mparticles_find_block_indices_ids_total(cuda_mprts);
  cuda_mparticles_dump(cuda_mprts);

  thrust::device_ptr<unsigned int> d_bidx(cuda_mprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cuda_mprts->d_id);
  thrust::stable_sort_by_key(d_bidx, d_bidx + cuda_mprts->nr_prts, d_id);
  cuda_mparticles_dump(cuda_mprts);

  cuda_mparticles_reorder_and_offsets(cuda_mprts);
  cuda_mparticles_dump(cuda_mprts);
  
  cuda_mparticles_free(cuda_mprts);
}
