
#include "cuda_mparticles.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

//#undef _GLIBCXX_USE_INT128

#include "cuda_sort2.h"
//#include "cuda_sort2_spine.h"
#include "particles_cuda.h"

#include <mrc_profile.h>

#if 0

#define STRIDE (9)
#define S_OOB (10)

#include <b40c/radixsort_reduction_kernel.h>
#include <b40c/radixsort_spine_kernel.h>
#include <b40c/radixsort_scanscatter_kernel.h>

using namespace b40c_thrust;

// OPT: dynamic shared mem

typedef unsigned int K;
typedef unsigned int V;

// blockIdx to rel offset
template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
struct PreShiftFunctor {
  __device__ __host__ __forceinline__ void operator()(K &converted_key) {
    converted_key = block_idx_to_dir1<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>(blockIdx.x, converted_key);
  }
  __device__ __host__ __forceinline__ static bool MustApply(){ return true;}
};

__global__ static void
set_new_offsets(int *offsets, int *bb_sums, int nr_blocks)
{
  int b = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  
  if (b <= nr_blocks) {
    offsets[b] = bb_sums[b * STRIDE + 0];
  }
}

static void
cuda_set_new_offsets(int *d_offsets, int *d_bb_sums, int nr_blocks)
{
  int dimGrid = (nr_blocks + 1 + THREADS_PER_BLOCK - 1) 
    / THREADS_PER_BLOCK;

  set_new_offsets<<<dimGrid, THREADS_PER_BLOCK>>>
    (d_offsets, d_bb_sums, nr_blocks);
}

class SortPairs {
  static const int RADIX_BITS = 4;
  static const int RADIX_DIGITS = 1 << RADIX_BITS;

  int _nr_blocks;
  thrust::device_vector<int> _bb_cnts;

  static int nr_blocks(const int b_mx[3])
  {
    return b_mx[0] * b_mx[1] * b_mx[2];
  }
  
public:
  int _b_mx[3];

  SortPairs(const int b_mx[3])
    : _nr_blocks(nr_blocks(b_mx)),
      _bb_cnts(9 * (_nr_blocks + 1))
  {
    _b_mx[0] = b_mx[0];
    _b_mx[1] = b_mx[1];
    _b_mx[2] = b_mx[2];

    int h_dirs[9][9][2];
    int h_dirs_inv[9][15];
    const int d_map[3][3] = {
      {  0, +1, -1, },
      { -1,  0, +1, },
      { +1, -1,  0, },
    };
    
    for (int d0 = 0; d0 < 9; d0++) {
      int d0y = d0 % 3 - 1;
      int d0z = d0 / 3 - 1;
      for (int e = 0; e < 9; e++) {
	int ey = d_map[d0 % 3][e % 3];
	int ez = d_map[d0 / 3][e / 3];
	
	int n = (1 - ez) << 2 | (1 - ey);

	if (d0y < 0 && ey < 0) ey = b_mx[1] - 1;
	if (d0y > 0 && ey > 0) ey = 1 - b_mx[1];
	if (d0z < 0 && ez < 0) ez = b_mx[2] - 1;
	if (d0z > 0 && ez > 0) ez = 1 - b_mx[2];

	h_dirs[d0][e][0] = n;
	h_dirs[d0][e][1] = ez * b_mx[1] + ey;
	h_dirs_inv[d0][n] = e;
      }
    }
    
    check(cudaMemcpyToSymbol(dirs, h_dirs, sizeof(dirs)));
    check(cudaMemcpyToSymbol(dirs_inv, h_dirs_inv, sizeof(dirs_inv)));
    // seems constant mem exists on the host side, too
    memcpy(dirs, h_dirs, sizeof(dirs));
    memcpy(dirs_inv, h_dirs_inv, sizeof(dirs_inv));
  }
  
  // ----------------------------------------------------------------------
  // sort

  template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
  void sort(unsigned int *d_bidx, unsigned int *d_alt_ids, int n, int *d_offsets)
  {
    static int pr_A, pr_B, pr_C;
    if (!pr_A) {
      pr_A = prof_register("sort_bottom_sum", 1., 0, 0);
      pr_B = prof_register("sort_top_scan", 1., 0, 0);
      pr_C = prof_register("sort_bottom_scan", 1., 0, 0);
    }
    
    prof_start(pr_A);
    reduction<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z> (d_offsets, d_bidx);
    prof_stop(pr_A);
    
    prof_start(pr_B);
    thrust::exclusive_scan(_bb_cnts.begin(), _bb_cnts.end(), _bb_cnts.begin());
    cuda_sync_if_enabled();
    prof_stop(pr_B);
    
    prof_start(pr_C);
    scan<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z> (d_offsets, d_bidx, d_alt_ids);
    cuda_set_new_offsets(d_offsets, thrust::raw_pointer_cast(&_bb_cnts[0]),
			 _nr_blocks);
    prof_stop(pr_C);
  }

private:
  // ----------------------------------------------------------------------
  // reduction

  template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
  void reduction(int *d_offsets, unsigned int *d_bidx)
  {
    const int threads = B40C_RADIXSORT_THREADS;
    RakingReduction2x<K, V, 0, RADIX_BITS, 0, 
		      PreShiftFunctor<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>,
		      NBLOCKS_Y, NBLOCKS_Z>
      <<<_nr_blocks, threads>>>
      (thrust::raw_pointer_cast(&_bb_cnts[0]), d_bidx, d_offsets);
    cuda_sync_if_enabled();
  }

  template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
  void scan(int *d_offsets, unsigned int *d_bidx, unsigned int *d_alt_ids)
  {
    const int threads = B40C_RADIXSORT_THREADS;

    ScanScatterDigits2x<K, V, 0, RADIX_BITS, 0,
			PreShiftFunctor<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>,
			NopFunctor<K>,
			NBLOCKS_Y, NBLOCKS_Z>
      <<<_nr_blocks, threads>>>
      (thrust::raw_pointer_cast(&_bb_cnts[0]), d_bidx, d_alt_ids, d_offsets);
    cuda_sync_if_enabled();
  }
};

// ======================================================================
// C interface

EXTERN_C void *
sort_pairs_create(const int b_mx[3])
{
  SortPairs *sp = new SortPairs(b_mx);
  
  return (void *) sp;
}

EXTERN_C void
sort_pairs_destroy(void *_sp)
{
  SortPairs *sp = (SortPairs *) _sp;

  delete sp;
}

EXTERN_C void
sort_pairs_device_2(void *_sp, unsigned int *d_bidx, unsigned int *d_alt_ids,
		    int n_part, int *d_offsets)
{
  SortPairs *sp = (SortPairs *) _sp;
  int *b_mx = sp->_b_mx;

#ifdef FAST_COMPILE
  if (b_mx[0] == _NBLOCKS_X && b_mx[1] == _NBLOCKS_Y && b_mx[2] == _NBLOCKS_Z) {
    sp->sort<_NBLOCKS_X, _NBLOCKS_Y, _NBLOCKS_Z>
      (d_bidx, d_alt_ids, n_part, d_offsets);
  } else {
    fprintf(stderr, "n_blocks %d x %d x %d not handled!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
#else
  if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    sp->sort<1, 8, 8> (d_bidx, d_alt_ids, n_part, d_offsets);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    sp->sort<1, 16, 16> (d_bidx, d_alt_ids, n_part, d_offsets);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    sp->sort<1, 32, 32> (d_bidx, d_alt_ids, n_part, d_offsets);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    sp->sort<1, 64, 64> (d_bidx, d_alt_ids, n_part, d_offsets);
  } else {
    fprintf(stderr, "n_blocks %d x %d x %d not handled!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
#endif
}

#endif

