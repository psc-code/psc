
#include "cuda_mfields.h"
#include "psc_cuda.h"

// the loops include 2 levels of ghost cells
// they really only need -1:2 and -1:1, respectively (for 1st order)
// but always doing 2:2 seems cheap enough

//#define BND 2

// FIXME, merge with F3_DEV{,_YZ}, OPT (precalc offset)

#define X3_DEV_OFF_YZ(fldnr, jy,jz)					\
  ((((fldnr)								\
     *mz + ((jz)+2))							\
    *my + ((jy)+2))							\
   *1 + (0))

#undef F3_DEV

#define F3_DEV(fldnr,ix,jy,jz)			\
  (d_flds)[X3_DEV_OFF_YZ(fldnr, jy,jz)]

#define F3_DDEV(d_flds, fldnr,ix,jy,jz)		\
  (d_flds)[X3_DEV_OFF_YZ(fldnr, jy,jz)]

__global__ static void
push_fields_E_yz(real *d_flds0, real dt, real cny, real cnz, int my, int mz,
		 unsigned int size, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < my - 2 * (2-BND) && iz < mz - 2 * (2-BND)))
    return;
  iy -= BND;
  iz -= BND;

  real *d_flds = d_flds0 + p * size;

  F3_DEV(EX, 0,iy,iz) +=
    cny * (F3_DEV(HZ, 0,iy,iz) - F3_DEV(HZ, 0,iy-1,iz)) -
    cnz * (F3_DEV(HY, 0,iy,iz) - F3_DEV(HY, 0,iy,iz-1)) -
    .5f * dt * F3_DEV(JXI, 0,iy,iz);
  
  F3_DEV(EY, 0,iy,iz) +=
    cnz * (F3_DEV(HX, 0,iy,iz) - F3_DEV(HX, 0,iy,iz-1)) -
    0.f -
    .5f * dt * F3_DEV(JYI, 0,iy,iz);
  
  F3_DEV(EZ, 0,iy,iz) +=
    0.f -
    cny * (F3_DEV(HX, 0,iy,iz) - F3_DEV(HX, 0,iy-1,iz)) -
    .5f * dt * F3_DEV(JZI, 0,iy,iz);
}

__global__ static void
push_fields_H_yz(real *d_flds0, real cny, real cnz, int my, int mz,
		 unsigned int size, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < my - 2 * (2-BND) && iz < mz - 2 * (2-BND)))
    return;
  iy -= BND;
  iz -= BND;

  real *d_flds = d_flds0 + p * size;

  F3_DEV(HX, 0,iy,iz) -=
    cny * (F3_DEV(EZ, 0,iy+1,iz) - F3_DEV(EZ, 0,iy,iz)) -
    cnz * (F3_DEV(EY, 0,iy,iz+1) - F3_DEV(EY, 0,iy,iz));
  
  F3_DEV(HY, 0,iy,iz) -=
    cnz * (F3_DEV(EX, 0,iy,iz+1) - F3_DEV(EX, 0,iy,iz)) -
    0.f;
  
  F3_DEV(HZ, 0,iy,iz) -=
    0.f -
    cny * (F3_DEV(EX, 0,iy+1,iz) - F3_DEV(EX, 0,iy,iz));
}

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

EXTERN_C void
cuda_push_fields_E_yz(struct cuda_mfields *cmflds, real dt)
{
  if (cmflds->n_patches == 0) {
    return;
  }

  real cny = .5f * dt / cmflds->dx[1];
  real cnz = .5f * dt / cmflds->dx[2];
  assert(cmflds->ldims[0] == 1);

  unsigned int size = cmflds->n_fields * cmflds->n_cells_per_patch;
  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (cmflds->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] * cmflds->n_patches };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_E_yz, (cmflds->d_flds, dt, cny, cnz, my, mz,
				size, grid[1]));
}

EXTERN_C void
cuda_push_fields_H_yz(struct cuda_mfields *cmflds, real dt)
{
  if (cmflds->n_patches == 0) {
    return;
  }

  real cny = .5f * dt / cmflds->dx[1];
  real cnz = .5f * dt / cmflds->dx[2];

  unsigned int size = cmflds->n_fields * cmflds->n_cells_per_patch;
  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (cmflds->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] * cmflds->n_patches };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_H_yz, (cmflds->d_flds, cny, cnz, my, mz,
				size, grid[1]));
}

EXTERN_C void
cuda_marder_correct_yz_gold(struct psc_mfields *mflds, struct psc_mfields *mf,
			    int p, int ldims[3], float fac[3],
			    int ly[3], int ry[3],
			    int lz[3], int rz[3])
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;
  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  struct cuda_mfields *cmf = psc_mfields_cuda(mf)->cmflds;
  fields_single_t f = cuda_mfields_get_host_fields(cmf);
  
  cuda_mfields_copy_from_device(cmflds, p, flds, EX, EX + 3);
  cuda_mfields_copy_from_device(cmf, p, f, 0, 1);
  
  for (int iz = -1; iz < ldims[2]; iz++) {
    for (int iy = -1; iy < ldims[1]; iy++) {
      if (iy >= -ly[1] && iy < ry[1] &&
	  iz >= -ly[2] && iz < ry[2]) {
	_F3_S(flds, EY, 0,iy,iz) += 
	  fac[1] * (_F3_S(f, 0, 0,iy+1,iz) - _F3_S(f, 0, 0,iy,iz));
	}
      
      if (iy >= -lz[1] && iy < rz[1] &&
	  iz >= -lz[2] && iz < rz[2]) {
	_F3_S(flds, EZ, 0,iy,iz) += 
	  fac[2] * (_F3_S(f, 0, 0,iy,iz+1) - _F3_S(f, 0, 0,iy,iz));
      }
    }
  }
  
  cuda_mfields_copy_to_device(cmflds, p, flds, EX, EX + 3);
  
  fields_single_t_dtor(&flds);
  fields_single_t_dtor(&f);
}

__global__ static void
marder_correct_yz(real *d_flds, real *d_f, float facy, float facz,
		  int lyy, int lyz, int ryy, int ryz,
		  int lzy, int lzz, int rzy, int rzz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  iy -= BND;
  iz -= BND;

  if (iy >= -lyy && iy < ryy &&
      iz >= -lyz && iz < ryz) {
    F3_DDEV(d_flds, EY, 0,iy,iz) += 
      facy * (F3_DDEV(d_f, 0, 0,iy+1,iz) - F3_DDEV(d_f, 0, 0,iy,iz));
  }
  
  if (iy >= -lzy && iy < rzy &&
      iz >= -lzz && iz < rzz) {
    F3_DDEV(d_flds, EZ, 0,iy,iz) += 
      facz * (F3_DDEV(d_f, 0, 0,iy,iz+1) - F3_DDEV(d_f, 0, 0,iy,iz));
  }
}

EXTERN_C void
cuda_marder_correct_yz(struct psc_mfields *mflds, struct psc_mfields *mf,
		       int p, int ldims[3], float fac[3],
		       int ly[3], int ry[3],
		       int lz[3], int rz[3])
{
#if 0
  cuda_marder_correct_yz_gold(mflds, mf, p, ldims, fac, ly, ry, lz, rz);
  return;
#endif

  if (mflds->nr_patches == 0) {
    return;
  }

  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  struct cuda_mfields *cmflds = mflds_cuda->cmflds;
  struct psc_mfields_cuda *mf_cuda = psc_mfields_cuda(mf);
  struct cuda_mfields *cmf = mf_cuda->cmflds;
  struct psc_patch *patch = &ppsc->patch[p];

  unsigned int size = cmflds->im[0] * cmflds->im[1] * cmflds->im[2];
  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] };

  RUN_KERNEL(dimGrid, dimBlock,
	     marder_correct_yz, (cmflds->d_flds + p * size * mflds->nr_fields,
				 cmf->d_flds + p * size * mf->nr_fields, fac[1], fac[2],
				 ly[1], ly[2], ry[1], ry[2],
				 lz[1], lz[2], rz[1], rz[2], my, mz));
}

// ======================================================================

__global__ static void
axpy_comp_yz(real *y_flds, int ym, float a, real *x_flds, int xm,
	     int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= my || iz >= mz) {
    return;
  }

  iy -= BND;
  iz -= BND;

  F3_DDEV(y_flds, ym, 0,iy,iz) += a * F3_DDEV(x_flds, xm, 0,iy,iz);
}

void
cuda_axpy_comp_yz(struct psc_mfields *mflds_y, int ym, float a, struct psc_mfields *mflds_x, int xm, int p)
{
  struct cuda_mfields *cmflds_y = psc_mfields_cuda(mflds_y)->cmflds;
  struct cuda_mfields *cmflds_x = psc_mfields_cuda(mflds_x)->cmflds;
  struct psc_patch *patch = &ppsc->patch[0];

  int my = cmflds_y->im[1];
  int mz = cmflds_y->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] };

  RUN_KERNEL(dimGrid, dimBlock,
	     axpy_comp_yz, (cmflds_y->d_flds_by_patch[p], ym, a, cmflds_x->d_flds_by_patch[p], xm, my, mz));
}

// ======================================================================
// FIXME, this should be more easily doable by just a cudaMemset()

__global__ static void
zero_comp_yz(real *x_flds, int xm, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= my || iz >= mz) {
    return;
  }

  iy -= BND;
  iz -= BND;

  F3_DDEV(x_flds, xm, 0,iy,iz) = 0.f;
}

void
cuda_zero_comp_yz(struct psc_mfields *mflds_x, int xm, int p)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds_x)->cmflds;
  struct psc_patch *patch = &ppsc->patch[0];

  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] };

  RUN_KERNEL(dimGrid, dimBlock,
	     zero_comp_yz, (cmflds->d_flds_by_patch[p], xm, my, mz));
}

// ======================================================================

__global__ static void
calc_dive_yz(real *flds, real *f, float dy, float dz,
	     int ldimsy, int ldimsz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= ldimsy || iz >= ldimsz) {
    return;
  }

  F3_DDEV(f, 0, 0,iy,iz) = 
    ((F3_DDEV(flds, EY, 0,iy,iz) - F3_DDEV(flds, EY, 0,iy-1,iz)) / dy +
     (F3_DDEV(flds, EZ, 0,iy,iz) - F3_DDEV(flds, EZ, 0,iy,iz-1)) / dz);
}

EXTERN_C void
cuda_calc_dive_yz(struct psc_mfields *mflds, struct psc_mfields *mf, int p)
{
  float dy = ppsc->patch[p].dx[1];
  float dz = ppsc->patch[p].dx[2];

  struct cuda_mfields *cmflds_flds = psc_mfields_cuda(mflds)->cmflds;
  struct cuda_mfields *cmflds_f = psc_mfields_cuda(mf)->cmflds;
  int *ldims = ppsc->patch[0].ldims;

  int my = cmflds_flds->im[1];
  int mz = cmflds_flds->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (ldims[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (ldims[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] };

  RUN_KERNEL(dimGrid, dimBlock,
	     calc_dive_yz, (cmflds_flds->d_flds_by_patch[p], cmflds_f->d_flds_by_patch[p], dy, dz, ldims[1], ldims[2], my, mz));
}

