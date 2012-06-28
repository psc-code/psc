
#define DIM_Z 1
#define DIM_YZ 2

__device__ static inline void
blockIdx_to_blockCrd(int bidx, int bi[3])
{
  bi[2] = bidx / (d_consts.b_mx[1] * d_consts.b_mx[0]);
  bidx -= bi[2] * (d_consts.b_mx[1] * d_consts.b_mx[0]);
  bi[1] = bidx / d_consts.b_mx[0];
  bidx -= bi[1] * d_consts.b_mx[0];
  bi[0] = bidx;
}

__device__ static inline int
blockPos_to_blockIdx(int block_pos[3])
{
#if BLOCKSIZE_X == 1
#ifdef NO_CHECKERBOARD
  return block_pos[2] * d_consts.b_mx[1] + block_pos[1];
#else
  int dimy = d_consts.b_mx[1] >> 1;
  return
    ((block_pos[1] & 1) << 0) |
    ((block_pos[2] & 1) << 1) | 
    (((block_pos[2] >> 1) * dimy + (block_pos[1] >> 1)) << 2);
#endif
#else
#error TBD
#endif
}

__device__ static inline void
blockIdx_to_cellPos(particles_cuda_dev_t *d_particles, int bidx, int ci[3])
{
#if BLOCKSIZE_X == 1
  int dimy = d_consts.b_mx[1] >> 1;
  int block_pos[3];
  block_pos[0] = 0;
  block_pos[1] = bidx & 1;
  block_pos[2] = (bidx >> 1) & 1;
  block_pos[1] |= ((bidx >> 2) % dimy) << 1;
  block_pos[2] |= ((bidx >> 2) / dimy) << 1;
  ci[0] = block_pos[0] * BLOCKSIZE_X;
  ci[1] = block_pos[1] * BLOCKSIZE_Y;
  ci[2] = block_pos[2] * BLOCKSIZE_Z;
#else
#error TBD
#endif
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
    real pos = xi[d] * d_consts.dxi[d] + shift;
    j[d] = cuda_fint(pos);
  }
}

__device__ static void
find_idx(const real xi[3], int j[3], real shift)
{
  int d;
  for (d = 0; d < 3; d++) {
    real pos = xi[d] * d_consts.dxi[d] + shift;
    j[d] = cuda_nint(pos);
  }
}

__device__ static void
find_idx_off(const real xi[3], int j[3], real h[3], real shift)
{
  int d;
  for (d = 0; d < 3; d++) {
    real pos = xi[d] * d_consts.dxi[d] + shift;
    j[d] = cuda_nint(pos);
    h[d] = j[d] - pos;
  }
}

