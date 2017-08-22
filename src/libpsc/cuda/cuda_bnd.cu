
#include "psc_cuda.h"
#include "cuda_mfields.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define PFX(x) cuda_bnd_##x
#include "constants.c"

#define SW (2) // FIXME

// OPT lots of optimization opportunity in the single-proc/patch ones,
// but they may not be that important for real production

__global__ static void
fill_ghosts_periodic_yz(real *d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_consts.mx[1] && iz < d_consts.mx[2]))
    return;

  bool inside = true;
  int jy = iy, jz = iz;
  if (jy < SW                  ) { jy += d_consts.mx[1] - 2*SW; inside = false; }
  if (jy >= d_consts.mx[1] - SW) { jy -= d_consts.mx[1] - 2*SW; inside = false; }
  if (jz < SW                  ) { jz += d_consts.mx[2] - 2*SW; inside = false; }
  if (jz >= d_consts.mx[2] - SW) { jz -= d_consts.mx[2] - 2*SW; inside = false; }

  if (inside)
    return;

  for (int m = mb; m < me; m++) {
    F3_DEV(m, 0,iy-SW,iz-SW) = F3_DEV(m, 0,jy-SW,jz-SW);
  }
}

EXTERN_C void
cuda_fill_ghosts_periodic_yz(int p, struct psc_fields *pf, int mb, int me)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(pf->mflds)->cmflds;
  cuda_bnd_set_constants(NULL, pf);

  struct psc_patch *patch = &ppsc->patch[p];
  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int dimGrid[2]  = { (patch->ldims[1] + 2*SW + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		      (patch->ldims[2] + 2*SW + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  RUN_KERNEL(dimGrid, dimBlock,
	     fill_ghosts_periodic_yz, (cmflds->d_flds_by_patch[p], mb, me));
}

__global__ static void
fill_ghosts_periodic_z(real *d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_consts.mx[1] && iz < d_consts.mx[2]))
    return;

  bool inside = true;
  int jy = iy, jz = iz;
  if (jz < SW                  ) { jz += d_consts.mx[2] - 2*SW; inside = false; }
  if (jz >= d_consts.mx[2] - SW) { jz -= d_consts.mx[2] - 2*SW; inside = false; }

  if (inside)
    return;

  for (int m = mb; m < me; m++) {
    F3_DEV(m, 0,iy-SW,iz-SW) = F3_DEV(m, 0,jy-SW,jz-SW);
  }
}

EXTERN_C void
cuda_fill_ghosts_periodic_z(int p, struct psc_fields *pf, int mb, int me)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(pf->mflds)->cmflds;
  cuda_bnd_set_constants(NULL, pf);

  struct psc_patch *patch = &ppsc->patch[p];
  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int dimGrid[2]  = { (patch->ldims[1] + 2*SW + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		      (patch->ldims[2] + 2*SW + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  RUN_KERNEL(dimGrid, dimBlock,
	     fill_ghosts_periodic_z, (cmflds->d_flds_by_patch[p], mb, me));
}

__global__ static void
add_ghosts_periodic_yz(real *d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_consts.mx[1] - 2*SW && iz < d_consts.mx[2] - 2*SW))
    return;

  if (iy < SW) {
    int jy = iy + (d_consts.mx[1] - 2*SW);
    int jz = iz;
    for (int m = mb; m < me; m++) {
      F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
    }
    if (iz < SW) {
      jz = iz + (d_consts.mx[2] - 2*SW);
      for (int m = mb; m < me; m++) {
	F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
      }
    }
    if (iz >= d_consts.mx[2] - 3*SW) {
      jz = iz - (d_consts.mx[2] - 2*SW);
      for (int m = mb; m < me; m++) {
	F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
      }
    }
  }
  if (iy >= d_consts.mx[1] - 3*SW) {
    int jy = iy - (d_consts.mx[1] - 2*SW);
    int jz = iz;
    for (int m = mb; m < me; m++) {
      F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
    }
    if (iz < SW) {
      jz = iz + (d_consts.mx[2] - 2*SW);
      for (int m = mb; m < me; m++) {
	F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
      }
    }
    if (iz >= d_consts.mx[2] - 3*SW) {
      jz = iz - (d_consts.mx[2] - 2*SW);
      for (int m = mb; m < me; m++) {
	F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
      }
    }
  }
  if (iz < SW) {
    int jy = iy, jz = iz + (d_consts.mx[2] - 2*SW);
    for (int m = mb; m < me; m++) {
      F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
    }
  }
  if (iz >= d_consts.mx[2] - 3*SW) {
    int jy = iy, jz = iz - (d_consts.mx[2] - 2*SW);
    for (int m = mb; m < me; m++) {
      F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
    }
  }
}

EXTERN_C void
cuda_add_ghosts_periodic_yz(int p, struct psc_fields *pf, int mb, int me)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(pf->mflds)->cmflds;
  cuda_bnd_set_constants(NULL, pf);

  struct psc_patch *patch = &ppsc->patch[p];
  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int dimGrid[2]  = { (patch->ldims[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		      (patch->ldims[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  RUN_KERNEL(dimGrid, dimBlock,
	     add_ghosts_periodic_yz, (cmflds->d_flds_by_patch[p], mb, me));
}

__global__ static void
add_ghosts_periodic_z(real *d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_consts.mx[1] - 2*SW && iz < d_consts.mx[2] - 2*SW))
    return;

  if (iz < SW) {
    int jy = iy, jz = iz + (d_consts.mx[2] - 2*SW);
    for (int m = mb; m < me; m++) {
      F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
    }
  }
  if (iz >= d_consts.mx[2] - 3*SW) {
    int jy = iy, jz = iz - (d_consts.mx[2] - 2*SW);
    for (int m = mb; m < me; m++) {
      F3_DEV(m, 0,iy,iz) += F3_DEV(m, 0,jy,jz);
    }
  }
}

EXTERN_C void
cuda_add_ghosts_periodic_z(int p, struct psc_fields *pf, int mb, int me)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(pf->mflds)->cmflds;
  cuda_bnd_set_constants(NULL, pf);

  struct psc_patch *patch = &ppsc->patch[p];
  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int dimGrid[2]  = { (patch->ldims[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		      (patch->ldims[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  RUN_KERNEL(dimGrid, dimBlock,
	     add_ghosts_periodic_z, (cmflds->d_flds_by_patch[p], mb, me));
}

template<bool lo, bool hi>
__global__ static void
conducting_wall_H_y(real *d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  if (iz >= d_consts.mx[2] - SW)
    return;

  int my = d_consts.mx[1] - 2*SW;

  if (lo) {
    F3_DEV(HY, 0,-1,iz) =  F3_DEV(HY, 0, 1,iz);
    F3_DEV(HX, 0,-1,iz) = -F3_DEV(HX, 0, 0,iz);
    F3_DEV(HZ, 0,-1,iz) = -F3_DEV(HZ, 0, 0,iz);
  }

  if (hi) {
    F3_DEV(HY, 0,my+1,iz) =  F3_DEV(HY, 0,my-1,iz);
    F3_DEV(HX, 0,my  ,iz) = -F3_DEV(HX, 0,my-1,iz);
    F3_DEV(HZ, 0,my  ,iz) = -F3_DEV(HZ, 0,my-1,iz);
  }
}

template<bool lo, bool hi>
__global__ static void
conducting_wall_E_y(real *d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  if (iz >= d_consts.mx[2] - SW)
    return;

  int my = d_consts.mx[1] - 2*SW;

  if (lo) {
    F3_DEV(EX, 0, 0,iz) =  0.;
    F3_DEV(EX, 0,-1,iz) =  F3_DEV(EX, 0, 1,iz);
    F3_DEV(EY, 0,-1,iz) = -F3_DEV(EY, 0, 0,iz);
    F3_DEV(EZ, 0, 0,iz) =  0.;
    F3_DEV(EZ, 0,-1,iz) =  F3_DEV(EZ, 0, 1,iz);
  }

  if (hi) {
    F3_DEV(EX, 0,my  ,iz) = 0.;
    F3_DEV(EX, 0,my+1,iz) =  F3_DEV(EX, 0, my-1,iz);
    F3_DEV(EY, 0,my,iz)   = -F3_DEV(EY, 0, my-1,iz);
    F3_DEV(EZ, 0,my,iz) = 0.;
    F3_DEV(EZ, 0,my+1,iz) =  F3_DEV(EZ, 0, my-1,iz);
  }
}

template<bool lo, bool hi>
__global__ static void
conducting_wall_J_y(real *d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  if (iz >= d_consts.mx[2] - SW)
    return;

  int my = d_consts.mx[1] - 2*SW;

  if (lo) {
    F3_DEV(JYI, 0, 0,iz) -= F3_DEV(JYI, 0,-1,iz);
    F3_DEV(JYI, 0,-1,iz) = 0.;
    F3_DEV(JXI, 0, 1,iz) += F3_DEV(JXI, 0,-1,iz);
    F3_DEV(JXI, 0,-1,iz) = 0.;
    F3_DEV(JZI, 0, 1,iz) += F3_DEV(JZI, 0,-1,iz);
    F3_DEV(JZI, 0,-1,iz) = 0.;
  }

  if (hi) {
    F3_DEV(JYI, 0,my-1,iz) -= F3_DEV(JYI, 0,my,iz);
    F3_DEV(JYI, 0,my  ,iz) = 0.;
    F3_DEV(JXI, 0,my-1,iz) += F3_DEV(JXI, 0,my+1,iz);
    F3_DEV(JXI, 0,my+1,iz) = 0.;
    F3_DEV(JZI, 0,my-1,iz) += F3_DEV(JZI, 0,my+1,iz);
    F3_DEV(JZI, 0,my+1,iz) = 0.;
  }
}

template<bool lo, bool hi>
static void
cuda_conducting_wall_H_y(int p, struct psc_fields *pf)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(pf->mflds)->cmflds;
  cuda_bnd_set_constants(NULL, pf);

  int dimGrid  = (ppsc->patch[p].ldims[2] + 2*SW + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_H_y<lo, hi> <<<dimGrid, BLOCKSIZE_Z>>> (cmflds->d_flds_by_patch[p]);
  cuda_sync_if_enabled();			       
}

template<bool lo, bool hi>
static void
cuda_conducting_wall_E_y(int p, struct psc_fields *pf)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(pf->mflds)->cmflds;
  cuda_bnd_set_constants(NULL, pf);

  int dimGrid  = (ppsc->patch[p].ldims[2] + 2*SW + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_E_y<lo, hi> <<<dimGrid, BLOCKSIZE_Z>>> (cmflds->d_flds_by_patch[p]);
  cuda_sync_if_enabled();			       
}

template<bool lo, bool hi>
static void
cuda_conducting_wall_J_y(int p, struct psc_fields *pf)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(pf->mflds)->cmflds;
  cuda_bnd_set_constants(NULL, pf);

  int dimGrid  = (ppsc->patch[p].ldims[2] + 2*SW + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_J_y<lo, hi> <<<dimGrid, BLOCKSIZE_Z>>> (cmflds->d_flds_by_patch[p]);
  cuda_sync_if_enabled();			       
}

EXTERN_C void
cuda_conducting_wall_H_lo_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_H_y<true, false>(p, pf);
}

EXTERN_C void
cuda_conducting_wall_H_hi_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_H_y<false, true>(p, pf);
}

EXTERN_C void
cuda_conducting_wall_H_lo_hi_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_H_y<true, true>(p, pf);
}

EXTERN_C void
cuda_conducting_wall_E_lo_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_E_y<true, false>(p, pf);
}

EXTERN_C void
cuda_conducting_wall_E_hi_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_E_y<false, true>(p, pf);
}

EXTERN_C void
cuda_conducting_wall_E_lo_hi_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_E_y<true, true>(p, pf);
}

EXTERN_C void
cuda_conducting_wall_J_lo_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_J_y<true, false>(p, pf);
}

EXTERN_C void
cuda_conducting_wall_J_hi_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_J_y<false, true>(p, pf);
}

EXTERN_C void
cuda_conducting_wall_J_lo_hi_y(int p, struct psc_fields *pf)
{
  cuda_conducting_wall_J_y<true, true>(p, pf);
}

