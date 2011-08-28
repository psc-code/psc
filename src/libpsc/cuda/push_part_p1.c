
// ----------------------------------------------------------------------
// push_part_one
//
// push one particle

__device__ static void
push_part_one(int n, particles_cuda_dev_t d_particles, real *d_flds, int l0[3])
{
  struct d_particle p;
  LOAD_PARTICLE(p, d_particles, n);
  real vxi[3];

  // x^n, p^n -> x^(n+0.5), p^n
  
  calc_vxi(vxi, p);
  push_xi(&p, vxi, .5f * d_dt);

  // field interpolation

  int lh[3], lg[3];
  real oh[3], og[3];
  find_idx_off(p.xi, lh, oh, real(-.5));
  find_idx_off(p.xi, lg, og, real(0.));

  real exq, eyq, ezq, hxq, hyq, hzq;
  INTERPOLATE_FIELD(exq, EX, g, g);
  INTERPOLATE_FIELD(eyq, EY, h, g);
  INTERPOLATE_FIELD(ezq, EZ, g, h);
  INTERPOLATE_FIELD(hxq, HX, h, h);
  INTERPOLATE_FIELD(hyq, HY, g, h);
  INTERPOLATE_FIELD(hzq, HZ, h, g);

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  
  push_pxi_dt(&p, exq, eyq, ezq, hxq, hyq, hzq);

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 

  calc_vxi(vxi, p);
  push_xi(&p, vxi, .5f * d_dt);

  STORE_PARTICLE_POS(p, d_particles, n);
  STORE_PARTICLE_MOM(p, d_particles, n);
}


// ----------------------------------------------------------------------
// push_part_p1
//
// push particles

__global__ static void
push_part_p1(int n_particles, particles_cuda_dev_t d_part, real *d_flds)
{
  int tid = threadIdx.x, bid = blockIdx.x;
  int block_begin = d_part.offsets[bid];
  int block_end   = d_part.offsets[bid+1]; // FIXME for blocksize != 1
  int ci[3];

  // cache fields

  blockIdx_to_blockCrd(blockIdx.x, ci);
  ci[0] *= BLOCKSIZE_X;
  ci[1] *= BLOCKSIZE_Y;
  ci[2] *= BLOCKSIZE_Z;

  cache_fields(d_flds, ci);
  
  for (int n = block_begin + tid; n < block_end; n += THREADS_PER_BLOCK) {
    /* printf("bid %d tid %d %d:%d ci %d,%d,%d\n", */
    /* 	   bid, tid, block_begin, block_end, ci[0], ci[1], ci[2]); */
    push_part_one(n, d_part, d_flds, ci); // FIXME d_flds?
  }
}

