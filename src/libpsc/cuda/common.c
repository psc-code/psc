
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
blockIdx_to_cellPos(particles_cuda_dev_t *d_particles, int bidx, int ci[3])
{
  int cidx = bidx * (BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z);
  ci[0] = d_particles->c_pos[3*cidx + 0];
  ci[1] = d_particles->c_pos[3*cidx + 1];
  ci[2] = d_particles->c_pos[3*cidx + 2];
}

__device__ static inline void
cellIdx_to_cellCrd_rel(int cidx, int ci[3])
{
#if BLOCKSIZE_X == 1 && BLOCKSIZE_Y == 4 && BLOCKSIZE_Z == 4
  ci[0] = 0;
  ci[1] = ((cidx & 4) >> 1) | ((cidx & 1) >> 0);
  ci[2] = ((cidx & 8) >> 2) | ((cidx & 2) >> 1);
#elif BLOCKSIZE_X == 1 && BLOCKSIZE_Y == 1 && BLOCKSIZE_Z == 1
  ci[0] = 0;
  ci[1] = 0;
  ci[2] = 0;
#else
#error TBD
#endif
}

__device__ static void
find_idx_1st(const real xi[3], int j[3], real shift)
{
  int d;
  for (d = 0; d < 3; d++) {
    real pos = xi[d] * d_dxi[d] + shift;
    j[d] = cuda_fint(pos);
  }
}

__device__ static void
find_idx_off_1st(const real xi[3], int j[3], real h[3], real shift)
{
  int d;
  for (d = 0; d < 3; d++) {
    real pos = xi[d] * d_dxi[d] + shift;
    j[d] = cuda_fint(pos);
    h[d] = pos - j[d];
  }
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

__device__ static real
ip1_to_grid_0(real h)
{
  return real(1.) - h;
}

__device__ static real
ip1_to_grid_p(real h)
{
  return h;
}

