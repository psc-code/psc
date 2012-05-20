
#include "psc_cuda.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define DIM DIM_YZ
#define PFX(x) yz_b_ ## x

#include "constants.c"
#include "common.c"
#include "common_push.c"
#include "common_fld_cache.c"

// ======================================================================

#define ARR_3_off(arr, d, o) (arr[(d-1)*3 + ((o)+1)])

__device__ static void
find_interp_to_grid_coeff(real *arr, const real h[3])
{
  for (int d = 1; d < 3; d++) {
    ARR_3_off(arr, d, -1) = ip_to_grid_m(h[d]);
    ARR_3_off(arr, d,  0) = ip_to_grid_0(h[d]);
    ARR_3_off(arr, d, +1) = ip_to_grid_p(h[d]);
  }
}

__device__ static real
get_ip_coeff_b(const real *arr, int d, int off)
{
  return ARR_3_off(arr, d, off);
}

__device__ static real
get_ip_coeff_b3(int off, const real h[3], int d)
{
  if (off == -1) {
    return ip_to_grid_m(h[d]);
  } else if (off == 0) {
    return ip_to_grid_0(h[d]);
  } else {
    return ip_to_grid_p(h[d]);
  }
}

__device__ static void
push_part_yz_b_one(int n, particles_cuda_dev_t d_particles, real *d_flds)
{
  struct d_particle p;
  LOAD_PARTICLE(p, d_particles, n);
  real h[3], vxi[3];

  // x^n, p^n -> x^(n+0.5), p^n
  
  calc_vxi(vxi, p);
  push_xi(&p, vxi, .5f * d_consts.dt);

  real __hh[2*3];
  int l[3];
  find_idx_off(p.xi, l, h, real(-.5));
  find_interp_to_grid_coeff(__hh, h);

  int j[3];
  real __gg[2*3];
  find_idx_off(p.xi, j, h, real(0.));
  find_interp_to_grid_coeff(__gg, h);

  // field interpolation
  
#define INTERPOLATE_FIELD_B(exq, EX, gg1, gg2, j0, j1, j2)		\
  real exq = 0;								\
  for (int dz = -1; dz <= 1; dz++) {					\
    for (int dy = -1; dy <= 1; dy++) {					\
      exq += get_ip_coeff_b(__##gg1, 1, dy) * get_ip_coeff_b(__##gg2, 2, dz) * \
	F3_DEV(EX, j0[0], j1[1]+dy, j2[2]+dz);				\
    }									\
  } do {} while(0)

  INTERPOLATE_FIELD_B(exq, EX, gg, gg, j, j, j);
  INTERPOLATE_FIELD_B(eyq, EY, hh, gg, j, l, j);
  INTERPOLATE_FIELD_B(ezq, EZ, gg, hh, j, j, l);
  INTERPOLATE_FIELD_B(hxq, HX, hh, hh, j, l, l);
  INTERPOLATE_FIELD_B(hyq, HY, gg, hh, j, j, l);
  INTERPOLATE_FIELD_B(hzq, HZ, hh, gg, j, l, j);

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  
  push_pxi_dt(&p, exq, eyq, ezq, hxq, hyq, hzq);

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 

  calc_vxi(vxi, p);
  push_xi(&p, vxi, .5f * d_consts.dt);

  STORE_PARTICLE_POS(p, d_particles, n);
  STORE_PARTICLE_MOM(p, d_particles, n);
}

__device__ static void
push_part_yz_b3_one(int n, particles_cuda_dev_t d_particles, real *d_flds)
{
  struct d_particle p;
  LOAD_PARTICLE(p, d_particles, n);
  real h[3], g[3], vxi[3];

  // x^n, p^n -> x^(n+0.5), p^n
  
  calc_vxi(vxi, p);
  push_xi(&p, vxi, .5f * d_consts.dt);

  int l[3];
  find_idx_off(p.xi, l, h, real(-.5));

  int j[3];
  find_idx_off(p.xi, j, g, real(0.));

  // field interpolation
  
#define INTERPOLATE_FIELD_B3(exq, EX, g1, g2, j0, j1, j2)		\
  real exq = 0;								\
  for (int dz = -1; dz <= 1; dz++) {					\
    for (int dy = -1; dy <= 1; dy++) {					\
      exq += get_ip_coeff_b3(dy, g1, 1) * get_ip_coeff_b3(dz, g2, 2) *	\
	F3_DEV(EX, j0[0], j1[1]+dy, j2[2]+dz);				\
    }									\
  } do {} while(0)

  INTERPOLATE_FIELD_B3(exq, EX, g, g, j, j, j);
  INTERPOLATE_FIELD_B3(eyq, EY, h, g, j, l, j);
  INTERPOLATE_FIELD_B3(ezq, EZ, g, h, j, j, l);
  INTERPOLATE_FIELD_B3(hxq, HX, h, h, j, l, l);
  INTERPOLATE_FIELD_B3(hyq, HY, g, h, j, j, l);
  INTERPOLATE_FIELD_B3(hzq, HZ, h, g, j, l, j);

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  
  push_pxi_dt(&p, exq, eyq, ezq, hxq, hyq, hzq);

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 

  calc_vxi(vxi, p);
  push_xi(&p, vxi, .5f * d_consts.dt);

  STORE_PARTICLE_POS(p, d_particles, n);
  STORE_PARTICLE_MOM(p, d_particles, n);
}

__global__ static void
push_part_yz_b(int n_part, particles_cuda_dev_t d_part, float *d_flds, int stride)
{
  int n = threadIdx.x + blockDim.x * blockIdx.x;

  while (n < n_part) {
    push_part_yz_b_one(n, d_part, d_flds);
    n += stride;
  }
}

__global__ static void
push_part_yz_b3(int n_part, particles_cuda_dev_t d_part, float *d_flds, int stride)
{
  int n = threadIdx.x + blockDim.x * blockIdx.x;

  while (n < n_part) {
    push_part_yz_b3_one(n, d_part, d_flds);
    n += stride;
  }
}

EXTERN_C void
__cuda_push_part_yz_b(struct psc_particles *prts, fields_cuda_t *pf)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  const int threadsPerBlock = 128;
  const int gridSize = 256;
  int dimBlock[2]  = { threadsPerBlock, 1 };
  int dimGrid[2] = { gridSize, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_yz_b, (prts->n_part, cuda->d_part, pf->d_flds,
			      gridSize * threadsPerBlock));
}

EXTERN_C void
__cuda_push_part_yz_b3(struct psc_particles *prts, fields_cuda_t *pf)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  const int threadsPerBlock = 128;
  const int gridSize = 256;
  int dimBlock[2]  = { threadsPerBlock, 1 };
  int dimGrid[2] = { gridSize, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_yz_b3, (prts->n_part, cuda->d_part, pf->d_flds,
			       gridSize * threadsPerBlock));
}

