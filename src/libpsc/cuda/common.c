
#define DIM_Z 1
#define DIM_YZ 2

__device__ static inline void
blockIdx_to_blockCrd(int bidx, int bi[3])
{
  bi[2] = bidx / (d_b_mx[1] * d_b_mx[0]);
  bidx -= bi[2] * (d_b_mx[1] * d_b_mx[0]);
  bi[1] = bidx / d_b_mx[0];
  bidx -= bi[1] * d_b_mx[0];
  bi[0] = bidx;
}

__device__ static inline void
cellIdx_to_cellCrd(int cidx, int ci[3])
{
#if BLOCKSIZE_Y == 4 && BLOCKSIZE_Z == 4
  blockIdx_to_blockCrd(cidx >> 4, ci);
  ci[1] = (ci[1] << 2) | ((cidx & 4) >> 1) | (cidx & 1);
  ci[2] = (ci[2] << 2) | ((cidx & 8) >> 2) | ((cidx & 2) >> 1);
#else
#error TBD
#endif
}

__device__ static void
find_idx(const real xi[3], int j[3], real shift)
{
  int d;
  for (d = 0; d < 3; d++) {
    real pos = xi[d] * d_dxi[d] + shift;
    j[d] = cuda_nint(pos);
  }
}

__device__ static void
find_idx_off(const real xi[3], int j[3], real h[3], real shift)
{
  int d;
  for (d = 0; d < 3; d++) {
    real pos = xi[d] * d_dxi[d] + shift;
    j[d] = cuda_nint(pos);
    h[d] = j[d] - pos;
  }
}

__device__ static real
ip_to_grid_m(real h)
{
  return real(.5) * sqr(real(.5) + h);
}

__device__ static real
ip_to_grid_0(real h)
{
  return real(.75) - sqr(h);
}

__device__ static real
ip_to_grid_p(real h)
{
  return real(.5) * sqr(real(.5) - h);
}

