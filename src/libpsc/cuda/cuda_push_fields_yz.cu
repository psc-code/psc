
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

// FIXME, duplicated
// ----------------------------------------------------------------------
// macros to access C (host) versions of the fields

#define F3_OFF_CUDA(pf, fldnr, jx,jy,jz)				\
  ((((((fldnr)								\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_CUDA(h_flds, pf, fldnr, jx,jy,jz)	\
  (h_flds[F3_OFF_CUDA(pf, fldnr, jx,jy,jz)])

#else

#define F3_CUDA(pf, fldnr, jx,jy,jz)				\
  (*({int off = F3_OFF_CUDA(pf, fldnr, jx,jy,jz);			\
      assert(fldnr >= 0 && fldnr < (pf)->nr_comp);			\
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &(h_flds[off]);						\
    }))

#endif


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
cuda_push_fields_E_yz(struct psc_mfields *mflds)
{
  if (mflds->nr_patches == 0) {
    return;
  }

  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  struct psc_patch *patch = &ppsc->patch[0];

  real dt = ppsc->dt;
  real cny = .5f * ppsc->dt / patch->dx[1];
  real cnz = .5f * ppsc->dt / patch->dx[2];
  assert(patch->ldims[0] == 1);

  unsigned int size = mflds->nr_fields *
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];
  int my = mflds_cuda->im[1];
  int mz = mflds_cuda->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] * mflds->nr_patches };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_E_yz, (mflds_cuda->d_flds, dt, cny, cnz, my, mz,
				size, grid[1]));
}

EXTERN_C void
cuda_push_fields_H_yz(struct psc_mfields *mflds)
{
  if (mflds->nr_patches == 0) {
    return;
  }

  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  struct psc_patch *patch = &ppsc->patch[0];

  real cny = .5f * ppsc->dt / patch->dx[1];
  real cnz = .5f * ppsc->dt / patch->dx[2];

  unsigned int size = mflds->nr_fields *
    mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];
  int my = mflds_cuda->im[1];
  int mz = mflds_cuda->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] * mflds->nr_patches };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_H_yz, (mflds_cuda->d_flds, cny, cnz, my, mz,
				size, grid[1]));
}

EXTERN_C void
cuda_marder_correct_yz_gold(struct psc_mfields *mflds, struct psc_mfields *mf,
			    int p, int ldims[3], float fac[3],
			    int ly[3], int ry[3],
			    int lz[3], int rz[3])
{
  struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
  struct psc_fields *f = psc_mfields_get_patch(mf, p);

  float *hflds = (float *) malloc(flds->nr_comp * flds->im[0] * flds->im[1] * flds->im[2] * sizeof(*hflds));
  float *hf = (float *) malloc(f->nr_comp * flds->im[0] * flds->im[1] * flds->im[2] * sizeof(*hf));
  
  __fields_cuda_from_device(flds, hflds, EX, EX + 3);
  __fields_cuda_from_device(f, hf, 0, 1);
  
  for (int iz = -1; iz < ldims[2]; iz++) {
    for (int iy = -1; iy < ldims[1]; iy++) {
      if (iy >= -ly[1] && iy < ry[1] &&
	  iz >= -ly[2] && iz < ry[2]) {
	F3_CUDA(hflds, flds, EY, 0,iy,iz) += 
	  fac[1] * (F3_CUDA(hf, f, 0, 0,iy+1,iz) - F3_CUDA(hf, f, 0, 0,iy,iz));
	}
      
      if (iy >= -lz[1] && iy < rz[1] &&
	  iz >= -lz[2] && iz < rz[2]) {
	F3_CUDA(hflds, flds, EZ, 0,iy,iz) += 
	  fac[2] * (F3_CUDA(hf, f, 0, 0,iy,iz+1) - F3_CUDA(hf, f, 0, 0,iy,iz));
      }
    }
  }
  
  __fields_cuda_to_device(flds, hflds, EX, EX + 3);
  
  free(hflds);
  free(hf);
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
  struct psc_mfields_cuda *mf_cuda = psc_mfields_cuda(mf);
  struct psc_patch *patch = &ppsc->patch[p];

  unsigned int size = mflds_cuda->im[0] * mflds_cuda->im[1] * mflds_cuda->im[2];
  int my = mflds_cuda->im[1];
  int mz = mflds_cuda->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] };

  RUN_KERNEL(dimGrid, dimBlock,
	     marder_correct_yz, (mflds_cuda->d_flds + p * size * mflds->nr_fields,
				 mf_cuda->d_flds + p * size * mf->nr_fields, fac[1], fac[2],
				 ly[1], ly[2], ry[1], ry[2],
				 lz[1], lz[2], rz[1], rz[2], my, mz));
}
