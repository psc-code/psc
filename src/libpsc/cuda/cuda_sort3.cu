
#include "psc_cuda.h"
#include "cuda_sort2.h"

#include <thrust/scan.h>

__constant__ __device__ int dirs[9][9][2];
__constant__ __device__ int dirs_inv[9][15];

// 0 - 8 is the 3x3 block of surrounding cells, 9 has additional particles being added
#define S_NEW (9)
#define STRIDE (10)
#define S_OOB (10)

#include <b40c/radixsort_reduction_kernel.h>
#include <b40c/radixsort_spine_kernel.h>
#include <b40c/radixsort_scanscatter_kernel3.h>

using namespace b40c_thrust;

typedef unsigned int K;
typedef unsigned int V;

static const int RADIX_BITS = 4;
//static const int RADIX_DIGITS = 1 << RADIX_BITS;


// FIXME -> header
EXTERN_C void cuda_copy_bidx_to_dev(struct psc_particles *prts, unsigned int *d_bidx, unsigned int *h_bidx);

// blockIdx to rel offset
template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
struct xPreShiftFunctor {
  __device__ __host__ __forceinline__ void operator()(K &converted_key) {
    if (converted_key == NBLOCKS_X * NBLOCKS_Y * NBLOCKS_Z) {
      converted_key = S_OOB;
    } else {
      int b_diff = blockIdx.x - converted_key + NBLOCKS_Y + 1;
      int d1 = b_diff % NBLOCKS_Y;
      int d2 = b_diff / NBLOCKS_Y;
      converted_key = d2 * 3 + d1;
    }
  }
  __device__ __host__ __forceinline__ static bool MustApply(){ return true;}
};

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
__global__ static void
cuda_set_bn_cnts(int *d_bb_cnts, unsigned int *d_bn_cnts)
{
  unsigned int b = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (b >= NBLOCKS_X * NBLOCKS_Y * NBLOCKS_Z)
    return;

  d_bb_cnts[b * STRIDE + S_NEW] = d_bn_cnts[b];
}

__global__ static void
cuda_move_recvd(int *d_bb_sums,
		unsigned int *d_bidx, unsigned int *d_alt_bidx,
		unsigned int *d_alt_ids,
		int n_part, int n_part_prev)
{
  int i = n_part_prev + threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (i >= n_part)
    return;
  
  int idx = d_bidx[i] * STRIDE + S_NEW;
  d_alt_ids[d_bb_sums[idx] + d_alt_bidx[i]] = i;
}

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

class SortPairs3 {
  int _nr_blocks;
  thrust::device_vector<int> _bb_cnts;
  thrust::device_vector<int> _bb_sums;
    
  static int nr_blocks(const int b_mx[3])
  {
    return b_mx[0] * b_mx[1] * b_mx[2];
  }
  
public:
  int _b_mx[3];

  SortPairs3(const int b_mx[3]) :
    _nr_blocks(nr_blocks(b_mx)),
    _bb_cnts(_nr_blocks * (STRIDE+1)),
    _bb_sums(_nr_blocks * (STRIDE+1))
  {
    _b_mx[0] = b_mx[0];
    _b_mx[1] = b_mx[1];
    _b_mx[2] = b_mx[2];
  }

  // ----------------------------------------------------------------------
  // sort

  template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
  void sort(unsigned int *d_bidx, unsigned int *d_alt_bidx, unsigned int *d_alt_ids,
	    int n_part, int *d_offsets, int n_part_prev, unsigned int *bn_cnts)
  {
    static int pr_A, pr_B, pr_C;
    if (!pr_A) {
      pr_A = prof_register("sort_bottom_sum", 1., 0, 0);
      pr_B = prof_register("sort_top_scan", 1., 0, 0);
      pr_C = prof_register("sort_bottom_scan", 1., 0, 0);
    }
    
    prof_start(pr_A);
    reduction<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>
      (d_offsets, d_bidx, n_part, n_part_prev, bn_cnts);
    prof_stop(pr_A);
    
    prof_start(pr_B);
    thrust::exclusive_scan(_bb_cnts.begin(), _bb_cnts.end(), _bb_sums.begin());
    cuda_sync_if_enabled();
    prof_stop(pr_B);
    
    prof_start(pr_C);
    scan<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>
      (d_offsets, d_bidx, d_alt_bidx, d_alt_ids, n_part, n_part_prev);
    cuda_set_new_offsets(d_offsets, thrust::raw_pointer_cast(&_bb_sums[0]),
			 _nr_blocks);
    prof_stop(pr_C);
  }

private:

  // ----------------------------------------------------------------------
  // reduction

  template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
  void reduction(int *d_offsets, unsigned int *d_bidx,
		 int n_part, int n_part_prev, unsigned int *bn_cnts)
  {
    // OPT, mostly unneeded, could just cache
    thrust::fill(_bb_cnts.begin(), _bb_cnts.end(), 0);

    // OPT, could leave out real interior counts (== 0)
    // copy to device
    thrust::device_vector<unsigned int> d_bn_cnts(bn_cnts, bn_cnts + _nr_blocks);

    int dimGrid = (n_part + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    cuda_set_bn_cnts<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>
      <<<dimGrid, THREADS_PER_BLOCK>>> (thrust::raw_pointer_cast(&_bb_cnts[0]),
					thrust::raw_pointer_cast(&d_bn_cnts[0]));
    cuda_sync_if_enabled();
    
    // set rest of array by counting 
    const int threads = B40C_RADIXSORT_THREADS;
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      xPreShiftFunctor<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>,
		      NBLOCKS_Y, NBLOCKS_Z>
      <<<_nr_blocks, threads>>> (thrust::raw_pointer_cast(&_bb_cnts[0]),
				 d_bidx, d_offsets);
    cuda_sync_if_enabled();
  }

  // ----------------------------------------------------------------------
  // reduction_host

  template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
  void reduction_host(struct psc_particles *prts, unsigned int *d_bidx,
		      int n_part, int n_part_prev)
  {
    int *d_bb_cnts = thrust::raw_pointer_cast(&_bb_cnts[0]);

    unsigned int *offsets = (unsigned int *) malloc((_nr_blocks + 1) * sizeof(*offsets));
    unsigned int *bidx = (unsigned int *) malloc(prts->n_part * sizeof(*bidx));
    int *bb_cnts = (int *) malloc(_nr_blocks * (STRIDE+1) * sizeof(*bb_cnts));
    memset(bb_cnts, 0, _nr_blocks * (STRIDE+1) * sizeof(*bb_cnts));
    cudaMemset(d_bb_cnts, 0, _nr_blocks * (STRIDE+1) * sizeof(*bb_cnts));
    cuda_copy_offsets_from_dev(prts, offsets);
    cuda_copy_bidx_from_dev(prts, bidx, d_bidx);
    
    // go by block for the old particles
    for (unsigned int bb = 0; bb < _nr_blocks; bb++) {
      for (int i = offsets[bb]; i < offsets[bb+1]; i++) {
	int idx;
	if (bidx[i] == _nr_blocks) {
	  idx = _nr_blocks * STRIDE + bb;
	} else {
	  assert(bidx[i] < _nr_blocks);
	  int b_diff = bb - bidx[i] + NBLOCKS_Y + 1;
	  int d1 = b_diff % NBLOCKS_Y;
	  int d2 = b_diff / NBLOCKS_Y;
	  idx = bidx[i] * STRIDE + (d2 * 3 + d1);
	}
	bb_cnts[idx]++;
      }
    }
    assert(offsets[_nr_blocks] == n_part_prev);
    
    // then do the new ones
    for (int i = n_part_prev; i < prts->n_part; i++) {
      assert(bidx[i] < _nr_blocks);
      bb_cnts[bidx[i] * STRIDE + S_NEW]++;
    }
    
    check(cudaMemcpy(d_bb_cnts, bb_cnts, _nr_blocks * (STRIDE+1) * sizeof(*bb_cnts),
		     cudaMemcpyHostToDevice));
    free(offsets);
    free(bidx);
    free(bb_cnts);
  }

  template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
  void scan(int *d_offsets, unsigned int *d_bidx,
	    unsigned int *d_alt_bidx, unsigned int *d_alt_ids,
	    int n_part, int n_part_prev)
  {
    const int threads = B40C_RADIXSORT_THREADS;
    
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			xPreShiftFunctor<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>,
			NopFunctor<K>,
			NBLOCKS_Y, NBLOCKS_Z> 
    <<<_nr_blocks, threads>>>
      (thrust::raw_pointer_cast(&_bb_sums[0]), d_bidx, d_alt_ids, d_offsets);
    cuda_sync_if_enabled();
  
    int dimsGrid = (n_part - n_part_prev + THREADS_PER_BLOCK - 1) /
      THREADS_PER_BLOCK;
    cuda_move_recvd<<<dimsGrid, THREADS_PER_BLOCK>>>
      (thrust::raw_pointer_cast(&_bb_sums[0]),
       d_bidx, d_alt_bidx, d_alt_ids, n_part, n_part_prev);
    cuda_sync_if_enabled();
  }

};

// ======================================================================
// C interface

EXTERN_C void *
sort_pairs_3_create(const int b_mx[3])
{
  SortPairs3 *sp = new SortPairs3(b_mx);
  
  return (void *) sp;
}

EXTERN_C void
sort_pairs_3_destroy(void *_sp)
{
  SortPairs3 *sp = (SortPairs3 *) _sp;

  delete sp;
}

EXTERN_C void
sort_pairs_3_device(void *_sp, unsigned int *d_bidx,
		    unsigned int *d_alt_bidx, unsigned int *d_alt_ids,
		    int n_part, int *d_offsets,
		    int n_part_prev, unsigned int *bn_cnts)
{
  SortPairs3 *sp = (SortPairs3 *) _sp;
  int *b_mx = sp->_b_mx;
  
  if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    sp->sort<1, 8, 8>(d_bidx, d_alt_bidx, d_alt_ids,
		      n_part, d_offsets, n_part_prev, bn_cnts);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    sp->sort<1, 16, 16>(d_bidx, d_alt_bidx, d_alt_ids,
			n_part, d_offsets, n_part_prev, bn_cnts);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    sp->sort<1, 32, 32>(d_bidx, d_alt_bidx, d_alt_ids,
			n_part, d_offsets, n_part_prev, bn_cnts);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    sp->sort<1, 64, 64>(d_bidx, d_alt_bidx, d_alt_ids,
			n_part, d_offsets, n_part_prev, bn_cnts);
  } else if (b_mx[0] == 1 && b_mx[1] == 128 && b_mx[2] == 128) {
    sp->sort<1, 128, 128>(d_bidx, d_alt_bidx, d_alt_ids,
			  n_part, d_offsets, n_part_prev, bn_cnts);
  } else {
    printf("need to add support for _b_mx %d %d\n", b_mx[1], b_mx[2]);
    assert(0);
  }
}
