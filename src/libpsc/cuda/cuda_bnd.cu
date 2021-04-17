
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include "psc.h"
#include "fields.hxx"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define SW (2) // FIXME

template <bool lo, bool hi>
__global__ static void conducting_wall_H_y(DFields d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  auto _d_flds = make_Fields3d<dim_xyz>(d_flds);

  if (iz >= d_flds.storage().shape(2) - SW)
    return;

  int my = d_flds.storage().shape(1) - 2 * SW;

  if (lo) {
    _d_flds(HY, 0, -1, iz) = _d_flds(HY, 0, 1, iz);
    _d_flds(HX, 0, -1, iz) = -_d_flds(HX, 0, 0, iz);
    _d_flds(HZ, 0, -1, iz) = -_d_flds(HZ, 0, 0, iz);
  }

  if (hi) {
    _d_flds(HY, 0, my + 1, iz) = _d_flds(HY, 0, my - 1, iz);
    _d_flds(HX, 0, my, iz) = -_d_flds(HX, 0, my - 1, iz);
    _d_flds(HZ, 0, my, iz) = -_d_flds(HZ, 0, my - 1, iz);
  }
}

template <bool lo, bool hi>
__global__ static void conducting_wall_E_y(DFields d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  auto _d_flds = make_Fields3d<dim_xyz>(d_flds);

  if (iz >= d_flds.storage().shape(2) - SW)
    return;

  int my = d_flds.storage().shape(1) - 2 * SW;

  if (lo) {
    _d_flds(EX, 0, 0, iz) = 0.;
    _d_flds(EX, 0, -1, iz) = _d_flds(EX, 0, 1, iz);
    _d_flds(EY, 0, -1, iz) = -_d_flds(EY, 0, 0, iz);
    _d_flds(EZ, 0, 0, iz) = 0.;
    _d_flds(EZ, 0, -1, iz) = _d_flds(EZ, 0, 1, iz);
  }

  if (hi) {
    _d_flds(EX, 0, my, iz) = 0.;
    _d_flds(EX, 0, my + 1, iz) = _d_flds(EX, 0, my - 1, iz);
    _d_flds(EY, 0, my, iz) = -_d_flds(EY, 0, my - 1, iz);
    _d_flds(EZ, 0, my, iz) = 0.;
    _d_flds(EZ, 0, my + 1, iz) = _d_flds(EZ, 0, my - 1, iz);
  }
}

template <bool lo, bool hi>
__global__ static void conducting_wall_J_y(DFields d_flds)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  auto _d_flds = make_Fields3d<dim_xyz>(d_flds);

  if (iz >= d_flds.storage().shape(2) - SW)
    return;

  int my = d_flds.storage().shape(1) - 2 * SW;

  if (lo) {
    _d_flds(JYI, 0, 0, iz) -= _d_flds(JYI, 0, -1, iz);
    _d_flds(JYI, 0, -1, iz) = 0.;
    _d_flds(JXI, 0, 1, iz) += _d_flds(JXI, 0, -1, iz);
    _d_flds(JXI, 0, -1, iz) = 0.;
    _d_flds(JZI, 0, 1, iz) += _d_flds(JZI, 0, -1, iz);
    _d_flds(JZI, 0, -1, iz) = 0.;
  }

  if (hi) {
    _d_flds(JYI, 0, my - 1, iz) -= _d_flds(JYI, 0, my, iz);
    _d_flds(JYI, 0, my, iz) = 0.;
    _d_flds(JXI, 0, my - 1, iz) += _d_flds(JXI, 0, my + 1, iz);
    _d_flds(JXI, 0, my + 1, iz) = 0.;
    _d_flds(JZI, 0, my - 1, iz) += _d_flds(JZI, 0, my + 1, iz);
    _d_flds(JZI, 0, my + 1, iz) = 0.;
  }
}

template <bool lo, bool hi>
static void cuda_conducting_wall_H_y(struct cuda_mfields* cmflds, int p)
{
  int dimGrid = (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_H_y<lo, hi><<<dimGrid, BLOCKSIZE_Z>>>((*cmflds)[p]);
  cuda_sync_if_enabled();
}

template <bool lo, bool hi>
static void cuda_conducting_wall_E_y(struct cuda_mfields* cmflds, int p)
{
  int dimGrid = (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_E_y<lo, hi><<<dimGrid, BLOCKSIZE_Z>>>((*cmflds)[p]);
  cuda_sync_if_enabled();
}

template <bool lo, bool hi>
static void cuda_conducting_wall_J_y(struct cuda_mfields* cmflds, int p)
{
  int dimGrid = (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_J_y<lo, hi><<<dimGrid, BLOCKSIZE_Z>>>((*cmflds)[p]);
  cuda_sync_if_enabled();
}

void cuda_conducting_wall_H_lo_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_H_y<true, false>(cmflds, p);
}

void cuda_conducting_wall_H_hi_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_H_y<false, true>(cmflds, p);
}

void cuda_conducting_wall_H_lo_hi_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_H_y<true, true>(cmflds, p);
}

void cuda_conducting_wall_E_lo_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_E_y<true, false>(cmflds, p);
}

void cuda_conducting_wall_E_hi_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_E_y<false, true>(cmflds, p);
}

void cuda_conducting_wall_E_lo_hi_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_E_y<true, true>(cmflds, p);
}

void cuda_conducting_wall_J_lo_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_J_y<true, false>(cmflds, p);
}

void cuda_conducting_wall_J_hi_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_J_y<false, true>(cmflds, p);
}

void cuda_conducting_wall_J_lo_hi_y(struct cuda_mfields* cmflds, int p)
{
  cuda_conducting_wall_J_y<true, true>(cmflds, p);
}
