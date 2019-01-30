
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include "psc.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define SW (2) // FIXME

// OPT lots of optimization opportunity in the single-proc/patch ones,
// but they may not be that important for float production

__global__ static void
fill_ghosts_periodic_yz(DFields d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_flds.im(1) && iz < d_flds.im(2)))
    return;

  bool inside = true;
  int jy = iy, jz = iz;
  if (jy < SW                ) { jy += d_flds.im(1) - 2*SW; inside = false; }
  if (jy >= d_flds.im(1) - SW) { jy -= d_flds.im(1) - 2*SW; inside = false; }
  if (jz < SW                ) { jz += d_flds.im(2) - 2*SW; inside = false; }
  if (jz >= d_flds.im(2) - SW) { jz -= d_flds.im(2) - 2*SW; inside = false; }

  if (inside)
    return;

  for (int m = mb; m < me; m++) {
    d_flds(m, 0,iy-SW,iz-SW) = d_flds(m, 0,jy-SW,jz-SW);
  }
}

void
cuda_fill_ghosts_periodic_yz(struct cuda_mfields *cmflds, int p, int mb, int me)
{
  assert(cmflds->ib[1] == -SW);
  assert(cmflds->ib[2] == -SW);

  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid((cmflds->im[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  fill_ghosts_periodic_yz<<<dimGrid, dimBlock>>>((*cmflds)[p], mb, me);
  cuda_sync_if_enabled();
}

__global__ static void
fill_ghosts_periodic_z(DFields d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_flds.im(1) && iz < d_flds.im(2)))
    return;

  bool inside = true;
  int jy = iy, jz = iz;
  if (jz < SW                ) { jz += d_flds.im(2) - 2*SW; inside = false; }
  if (jz >= d_flds.im(2) - SW) { jz -= d_flds.im(2) - 2*SW; inside = false; }

  if (inside)
    return;

  for (int m = mb; m < me; m++) {
    d_flds(m, 0,iy-SW,iz-SW) = d_flds(m, 0,jy-SW,jz-SW);
  }
}

void
cuda_fill_ghosts_periodic_z(struct cuda_mfields *cmflds, int p, int mb, int me)
{
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid((cmflds->im[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  fill_ghosts_periodic_z<<<dimGrid, dimBlock>>>((*cmflds)[p], mb, me);
  cuda_sync_if_enabled();
}

__global__ static void
add_ghosts_periodic_yz(DFields d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_flds.im(1) - 2*SW && iz < d_flds.im(2) - 2*SW))
    return;

  if (iy < SW) {
    int jy = iy + (d_flds.im(1) - 2*SW);
    int jz = iz;
    for (int m = mb; m < me; m++) {
      d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
    }
    if (iz < SW) {
      jz = iz + (d_flds.im(2) - 2*SW);
      for (int m = mb; m < me; m++) {
	d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
      }
    }
    if (iz >= d_flds.im(2) - 3*SW) {
      jz = iz - (d_flds.im(2) - 2*SW);
      for (int m = mb; m < me; m++) {
	d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
      }
    }
  }
  if (iy >= d_flds.im(1) - 3*SW) {
    int jy = iy - (d_flds.im(1) - 2*SW);
    int jz = iz;
    for (int m = mb; m < me; m++) {
      d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
    }
    if (iz < SW) {
      jz = iz + (d_flds.im(2) - 2*SW);
      for (int m = mb; m < me; m++) {
	d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
      }
    }
    if (iz >= d_flds.im(2) - 3*SW) {
      jz = iz - (d_flds.im(2) - 2*SW);
      for (int m = mb; m < me; m++) {
	d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
      }
    }
  }
  if (iz < SW) {
    int jy = iy, jz = iz + (d_flds.im(2) - 2*SW);
    for (int m = mb; m < me; m++) {
      d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
    }
  }
  if (iz >= d_flds.im(2) - 3*SW) {
    int jy = iy, jz = iz - (d_flds.im(2) - 2*SW);
    for (int m = mb; m < me; m++) {
      d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
    }
  }
}

void
cuda_add_ghosts_periodic_yz(struct cuda_mfields *cmflds, int p, int mb, int me)
{
  Int3 ldims = cmflds->grid().ldims;
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid((ldims[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (ldims[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);

  add_ghosts_periodic_yz<<<dimGrid, dimBlock>>>((*cmflds)[p], mb, me);
  cuda_sync_if_enabled();
}

__global__ static void
add_ghosts_periodic_z(DFields d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_flds.im(1) - 2*SW && iz < d_flds.im(2) - 2*SW))
    return;

  if (iz < SW) {
    int jy = iy, jz = iz + (d_flds.im(2) - 2*SW);
    for (int m = mb; m < me; m++) {
      d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
    }
  }
  if (iz >= d_flds.im(2) - 3*SW) {
    int jy = iy, jz = iz - (d_flds.im(2) - 2*SW);
    for (int m = mb; m < me; m++) {
      d_flds(m, 0,iy,iz) += d_flds(m, 0,jy,jz);
    }
  }
}

void
cuda_add_ghosts_periodic_z(struct cuda_mfields *cmflds, int p, int mb, int me)
{
  Int3 ldims = cmflds->grid().ldims;
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid((ldims[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (ldims[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  add_ghosts_periodic_z<<<dimGrid, dimBlock>>>((*cmflds)[p], mb, me);
  cuda_sync_if_enabled();
}

template<bool lo, bool hi>
__global__ static void
conducting_wall_H_y(DFields d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  if (iz >= d_flds.im(2) - SW)
    return;

  int my = d_flds.im(1) - 2*SW;

  if (lo) {
    d_flds(HY, 0,-1,iz) =  d_flds(HY, 0, 1,iz);
    d_flds(HX, 0,-1,iz) = -d_flds(HX, 0, 0,iz);
    d_flds(HZ, 0,-1,iz) = -d_flds(HZ, 0, 0,iz);
  }

  if (hi) {
    d_flds(HY, 0,my+1,iz) =  d_flds(HY, 0,my-1,iz);
    d_flds(HX, 0,my  ,iz) = -d_flds(HX, 0,my-1,iz);
    d_flds(HZ, 0,my  ,iz) = -d_flds(HZ, 0,my-1,iz);
  }
}

template<bool lo, bool hi>
__global__ static void
conducting_wall_E_y(DFields d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  if (iz >= d_flds.im(2) - SW)
    return;

  int my = d_flds.im(1) - 2*SW;

  if (lo) {
    d_flds(EX, 0, 0,iz) =  0.;
    d_flds(EX, 0,-1,iz) =  d_flds(EX, 0, 1,iz);
    d_flds(EY, 0,-1,iz) = -d_flds(EY, 0, 0,iz);
    d_flds(EZ, 0, 0,iz) =  0.;
    d_flds(EZ, 0,-1,iz) =  d_flds(EZ, 0, 1,iz);
  }

  if (hi) {
    d_flds(EX, 0,my  ,iz) = 0.;
    d_flds(EX, 0,my+1,iz) =  d_flds(EX, 0, my-1,iz);
    d_flds(EY, 0,my,iz)   = -d_flds(EY, 0, my-1,iz);
    d_flds(EZ, 0,my,iz) = 0.;
    d_flds(EZ, 0,my+1,iz) =  d_flds(EZ, 0, my-1,iz);
  }
}

template<bool lo, bool hi>
__global__ static void
conducting_wall_J_y(DFields d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  if (iz >= d_flds.im(2) - SW)
    return;

  int my = d_flds.im(1) - 2*SW;

  if (lo) {
    d_flds(JYI, 0, 0,iz) -= d_flds(JYI, 0,-1,iz);
    d_flds(JYI, 0,-1,iz) = 0.;
    d_flds(JXI, 0, 1,iz) += d_flds(JXI, 0,-1,iz);
    d_flds(JXI, 0,-1,iz) = 0.;
    d_flds(JZI, 0, 1,iz) += d_flds(JZI, 0,-1,iz);
    d_flds(JZI, 0,-1,iz) = 0.;
  }

  if (hi) {
    d_flds(JYI, 0,my-1,iz) -= d_flds(JYI, 0,my,iz);
    d_flds(JYI, 0,my  ,iz) = 0.;
    d_flds(JXI, 0,my-1,iz) += d_flds(JXI, 0,my+1,iz);
    d_flds(JXI, 0,my+1,iz) = 0.;
    d_flds(JZI, 0,my-1,iz) += d_flds(JZI, 0,my+1,iz);
    d_flds(JZI, 0,my+1,iz) = 0.;
  }
}

template<bool lo, bool hi>
static void
cuda_conducting_wall_H_y(struct cuda_mfields *cmflds, int p)
{
  int dimGrid  = (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_H_y<lo, hi> <<<dimGrid, BLOCKSIZE_Z>>> ((*cmflds)[p]);
  cuda_sync_if_enabled();			       
}

template<bool lo, bool hi>
static void
cuda_conducting_wall_E_y(struct cuda_mfields *cmflds, int p)
{
  int dimGrid  = (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_E_y<lo, hi> <<<dimGrid, BLOCKSIZE_Z>>> ((*cmflds)[p]);
  cuda_sync_if_enabled();			       
}

template<bool lo, bool hi>
static void
cuda_conducting_wall_J_y(struct cuda_mfields *cmflds, int p)
{
  int dimGrid  = (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_J_y<lo, hi> <<<dimGrid, BLOCKSIZE_Z>>> ((*cmflds)[p]);
  cuda_sync_if_enabled();			       
}

void
cuda_conducting_wall_H_lo_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_H_y<true, false>(cmflds, p);
}

void
cuda_conducting_wall_H_hi_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_H_y<false, true>(cmflds, p);
}

void
cuda_conducting_wall_H_lo_hi_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_H_y<true, true>(cmflds, p);
}

void
cuda_conducting_wall_E_lo_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_E_y<true, false>(cmflds, p);
}

void
cuda_conducting_wall_E_hi_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_E_y<false, true>(cmflds, p);
}

void
cuda_conducting_wall_E_lo_hi_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_E_y<true, true>(cmflds, p);
}

void
cuda_conducting_wall_J_lo_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_J_y<true, false>(cmflds, p);
}

void
cuda_conducting_wall_J_hi_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_J_y<false, true>(cmflds, p);
}

void
cuda_conducting_wall_J_lo_hi_y(struct cuda_mfields *cmflds, int p)
{
  cuda_conducting_wall_J_y<true, true>(cmflds, p);
}

