

__device__ static void
reduce_sum_sdata()
{
  unsigned int tid = threadIdx.x;

  // adapted from CUDA SDK
  __syncthreads();
  if (THREADS_PER_BLOCK >= 512) {
    if (tid < 256) 
      forall_j(SDATA(tid,j) += SDATA(tid + 256,j););
    __syncthreads(); 
  }
  if (THREADS_PER_BLOCK >= 256) {
    if (tid < 128) 
      forall_j(SDATA(tid,j) += SDATA(tid + 128,j););
    __syncthreads(); 
  }
  if (THREADS_PER_BLOCK >= 128) {
    if (tid <  64)
      forall_j(SDATA(tid,j) += SDATA(tid + 64,j););
    __syncthreads();
  }
  
  if (tid < 32) {
    if (THREADS_PER_BLOCK >=  64) {
      forall_j(SDATA(tid,j) += SDATA(tid + 32,j););
    }
    if (THREADS_PER_BLOCK >=  32) {
      forall_j(SDATA(tid,j) += SDATA(tid + 16,j););
    }
    if (THREADS_PER_BLOCK >=  16) {
      forall_j(SDATA(tid,j) += SDATA(tid +  8,j););
    }
    if (THREADS_PER_BLOCK >=   8) {
      forall_j(SDATA(tid,j) += SDATA(tid +  4,j););
    }
    if (THREADS_PER_BLOCK >=   4) {
      forall_j(SDATA(tid,j) += SDATA(tid +  2,j););
    }
    if (THREADS_PER_BLOCK >=   2) {
      forall_j(SDATA(tid,j) += SDATA(tid +  1,j););
    }
  }
}

