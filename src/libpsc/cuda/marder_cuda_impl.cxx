
#include "marder_cuda_impl.hxx"

#define BND 2
#define BLOCKSIZE_X 8
#define BLOCKSIZE_Y 8
#define BLOCKSIZE_Z 8

void cuda_marder_correct_yz_gold(MfieldsStateCuda& mflds, MfieldsCuda& mf,
                                 int p, Float3 fac, Int3 ly, Int3 ry, Int3 lz,
                                 Int3 rz)
{
  auto h_mflds = hostMirror(mflds);
  auto h_mf = hostMirror(mf);

  copy(mflds, h_mflds);
  copy(mf, h_mf);

  auto flds = make_Fields3d<dim_xyz>(h_mflds[p]);
  auto f = make_Fields3d<dim_xyz>(h_mf[p]);

  Int3 ldims = h_mflds.grid().ldims;
  for (int iz = -1; iz < ldims[2]; iz++) {
    for (int iy = -1; iy < ldims[1]; iy++) {
      if (iy >= -ly[1] && iy < ry[1] && iz >= -ly[2] && iz < ry[2]) {
        flds(EY, 0, iy, iz) += fac[1] * (f(0, 0, iy + 1, iz) - f(0, 0, iy, iz));
      }

      if (iy >= -lz[1] && iy < rz[1] && iz >= -lz[2] && iz < rz[2]) {
        flds(EZ, 0, iy, iz) += fac[2] * (f(0, 0, iy, iz + 1) - f(0, 0, iy, iz));
      }
    }
  }

  copy(h_mflds, mflds);
}

void cuda_marder_correct_yz(MfieldsStateCuda& mflds, MfieldsCuda& mf, int p,
                            Float3 fac, Int3 ly, Int3 ry, Int3 lz, Int3 rz)
{
#if 0
  cuda_marder_correct_yz_gold(mflds, mf, p, fac, ly, ry, lz, rz);
  return;
#endif

  auto gt_flds = mflds.gt().view(_all, _all, _all, _all, p).to_kernel();
  auto gt_f = mf.gt().view(_all, _all, _all, 0, p).to_kernel();

  gt::launch<2>(
    {gt_flds.shape(1), gt_flds.shape(2)}, GT_LAMBDA(int iy, int iz) {
      if (iy - BND >= -ly[1] && iy - BND < ry[1] && iz - BND >= -ly[2] &&
          iz - BND < ry[2]) {
        gt_flds(0, iy, iz, EY) =
          gt_flds(0, iy, iz, EY) +
          fac[1] * (gt_f(0, iy + 1, iz) - gt_f(0, iy, iz));
      }

      if (iy - BND >= -lz[1] && iy - BND < rz[1] && iz - BND >= -lz[2] &&
          iz - BND < rz[2]) {
        gt_flds(0, iy, iz, EZ) =
          gt_flds(0, iy, iz, EZ) +
          fac[2] * (gt_f(0, iy, iz + 1) - gt_f(0, iy, iz));
      }
    });
  cuda_sync_if_enabled();
}

__global__ static void marder_correct_xyz(
  DFields d_flds, DFields d_f, float facx, float facy, float facz, int lxx,
  int lxy, int lxz, int rxx, int rxy, int rxz, int lyx, int lyy, int lyz,
  int ryx, int ryy, int ryz, int lzx, int lzy, int lzz, int rzx, int rzy,
  int rzz, int mx, int my, int mz)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  int iz = blockIdx.z * blockDim.z + threadIdx.z;

  ix -= BND;
  iy -= BND;
  iz -= BND;

  auto _d_flds = make_Fields3d<dim_xyz>(d_flds);
  auto _d_f = make_Fields3d<dim_xyz>(d_f);

  if (ix >= -lxx && ix < rxx && iy >= -lxy && iy < rxy && iz >= -lxz &&
      iz < rxz) {
    _d_flds(EX, ix, iy, iz) +=
      facx * (_d_f(0, ix + 1, iy, iz) - _d_f(0, ix, iy, iz));
  }

  if (ix >= -lyx && ix < ryx && iy >= -lyy && iy < ryy && iz >= -lyz &&
      iz < ryz) {
    _d_flds(EY, ix, iy, iz) +=
      facy * (_d_f(0, ix, iy + 1, iz) - _d_f(0, ix, iy, iz));
  }

  if (ix >= -lzx && ix < rzx && iy >= -lzy && iy < rzy && iz >= -lzz &&
      iz < rzz) {
    _d_flds(EZ, ix, iy, iz) +=
      facz * (_d_f(0, ix, iy, iz + 1) - _d_f(0, ix, iy, iz));
  }
}

void cuda_marder_correct_xyz(struct cuda_mfields* cmflds,
                             struct cuda_mfields* cmf, int p, float fac[3],
                             int lx[3], int rx[3], int ly[3], int ry[3],
                             int lz[3], int rz[3])
{
  if (cmflds->n_patches() == 0) {
    return;
  }

  int mx = cmflds->im(0);
  int my = cmflds->im(1);
  int mz = cmflds->im(2);

  int grid[3] = {(mx + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
                 (my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2]);

  marder_correct_xyz<<<dimGrid, dimBlock>>>(
    (*cmflds)[p], (*cmf)[p], fac[0], fac[1], fac[2], lx[0], lx[1], lx[2], rx[0],
    rx[1], rx[2], ly[0], ly[1], ly[2], ry[0], ry[1], ry[2], lz[0], lz[1], lz[2],
    rz[0], rz[1], rz[2], mx, my, mz);
  cuda_sync_if_enabled();
}
