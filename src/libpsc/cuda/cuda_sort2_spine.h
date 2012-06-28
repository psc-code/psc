
#include "cuda_sort2.h"

static __constant__ __device__ int dirs[9][9][2];
static __constant__ __device__ int dirs_inv[9][15];

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
static __host__ __device__ __forceinline__ void
spine_one(int b2, int *spine, int &sum, bool cnt, int d0)
{
  const int NBLOCKS = NBLOCKS_X * NBLOCKS_Y * NBLOCKS_Z;

  #pragma unroll
  for (int n = 0; n < 9; n++) {
    int d = dirs[d0][n][0];
    int b = b2 + dirs[d0][n][1];
    int val = spine[d * NBLOCKS + b];
    if (!cnt) spine[d * NBLOCKS + b] = sum;
    sum += val;
  }
}

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z, int ELEMENTS_PER_THREAD>
static void __host__ __device__
do_spine_tid(int *spine, int *sum, int tid, bool cnt)
{
  for (int b2 = tid * ELEMENTS_PER_THREAD; b2 < (tid + 1) * ELEMENTS_PER_THREAD; b2++) {
    int b2y = b2 % NBLOCKS_Y;
    int b2z = b2 / NBLOCKS_Y;
    int d0y, d0z;
    if (b2y == 0) {
      d0y = 0;
    } else if (b2y == NBLOCKS_Y - 1) {
      d0y = 2;
    } else {
      d0y = 1;
    }
    if (b2z == 0) {
      d0z = 0;
    } else if (b2z == NBLOCKS_Z - 1) {
      d0z = 2;
    } else {
      d0z = 1;
    }
    spine_one<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>(b2, spine, sum[tid], cnt, d0z * 3 + d0y);
  }
}

const int SPINE_NR_ELEMENT_THREADS = 128; // fermi: 512 ?

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
void
cuda_spine_host(thrust::device_vector<int> &d_spine)
{
  const int NR_ELEMENT_THREADS = SPINE_NR_ELEMENT_THREADS;
  const int NBLOCKS = NBLOCKS_X * NBLOCKS_Y * NBLOCKS_Z;
  const int ELEMENTS_PER_THREAD = NBLOCKS / NR_ELEMENT_THREADS;
  b40c_thrust::SuppressUnusedConstantWarning(ELEMENTS_PER_THREAD);
  assert(NBLOCKS % NR_ELEMENT_THREADS == 0);
  
  thrust::host_vector<int> h_spine = d_spine;

  static int pr;
  if (!pr) {
    pr = prof_register("sort_top_scan_c", 1., 0, 0);
  }
  
  prof_start(pr);
  int sum[NR_ELEMENT_THREADS];
  // reduce
  for (int tid = 0; tid < NR_ELEMENT_THREADS; tid++) {
    sum[tid] = 0;
    do_spine_tid<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z, ELEMENTS_PER_THREAD>(&h_spine[0], sum, tid, true);
  }

  // top-level
  int partial_sum = 0;
  for (int tid = 0; tid < NR_ELEMENT_THREADS; tid++) {
    int val = sum[tid];
    sum[tid] = partial_sum;
    partial_sum += val;
  }

  // scan
  for (int tid = 0; tid < NR_ELEMENT_THREADS; tid++) {
    do_spine_tid<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z, ELEMENTS_PER_THREAD>(&h_spine[0], sum, tid, false);
  }
  prof_stop(pr);
  
  d_spine = h_spine;
}

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
static __global__ void
cuda_spine(int *d_spine)
{
  const int NR_ELEMENT_THREADS = SPINE_NR_ELEMENT_THREADS;
  const int NBLOCKS = NBLOCKS_X * NBLOCKS_Y * NBLOCKS_Z;
  const int ELEMENTS_PER_THREAD = NBLOCKS / NR_ELEMENT_THREADS;
  b40c_thrust::SuppressUnusedConstantWarning(ELEMENTS_PER_THREAD);
  
  __shared__ int sum[NR_ELEMENT_THREADS];

  sum[threadIdx.x] = 0;
  do_spine_tid<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z, ELEMENTS_PER_THREAD>
    (d_spine, sum, threadIdx.x, true);

  __syncthreads();

  if (threadIdx.x == 0) {
    // top-level
    int partial_sum = 0;
    for (int tid = 0; tid < NR_ELEMENT_THREADS; tid++) {
      int val = sum[tid];
      sum[tid] = partial_sum;
      partial_sum += val;
    }
  }

  __syncthreads();

  do_spine_tid<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z, ELEMENTS_PER_THREAD>
    (d_spine, sum, threadIdx.x, false);
}

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
void
cuda_spine_device(thrust::device_vector<int> &d_spine)
{
  const int NBLOCKS = NBLOCKS_X * NBLOCKS_Y * NBLOCKS_Z;
  assert(NBLOCKS % SPINE_NR_ELEMENT_THREADS == 0);
  
  cuda_spine<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z><<<1, SPINE_NR_ELEMENT_THREADS>>>
    (thrust::raw_pointer_cast(&d_spine[0]));
  cuda_sync_if_enabled();

#if 0
  {
  printf("abc %d %d %d %d\n", dirs[0][0][0], dirs[0][0][1], dirs[0][1][0], dirs[0][1][1]);
  thrust::host_vector<int> h_spine = d_spine;
  for (int i = 0; i < 20; i++) {
    printf("%d: %d\n", i, h_spine[i]);
  }
  exit(0);
  }
#endif
}

