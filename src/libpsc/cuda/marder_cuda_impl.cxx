
#include "marder_cuda_impl.hxx"

#define BND 2
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

void cuda_marder_correct_yz_gold(struct cuda_mfields* cmflds,
                                 struct cuda_mfields* cmf, int p, float fac[3],
                                 int ly[3], int ry[3], int lz[3], int rz[3])
{
  auto mflds = hostMirror(*cmflds);
  auto mf = hostMirror(*cmf);

  copy(*cmflds, mflds);
  copy(*cmf, mf);

  auto flds = make_Fields3d<dim_xyz>(mflds[p]);
  auto f = make_Fields3d<dim_xyz>(mf[p]);

  Int3 ldims = cmflds->grid().ldims;
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

  copy(mflds, *cmflds);
}

template <typename E1, typename E2>
__global__ static void marder_correct_yz(E1 gt_flds, E2 gt_f, float facy,
                                         float facz, int lyy, int lyz, int ryy,
                                         int ryz, int lzy, int lzz, int rzy,
                                         int rzz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy - BND >= -lyy && iy - BND < ryy && iz - BND >= -lyz &&
      iz - BND < ryz) {
    gt_flds(0, iy, iz, EY) =
      gt_flds(0, iy, iz, EY) + facy * (gt_f(0, iy + 1, iz) - gt_f(0, iy, iz));
  }

  if (iy - BND >= -lzy && iy - BND < rzy && iz - BND >= -lzz &&
      iz - BND < rzz) {
    gt_flds(0, iy, iz, EZ) =
      gt_flds(0, iy, iz, EZ) + facz * (gt_f(0, iy, iz + 1) - gt_f(0, iy, iz));
  }
}

void cuda_marder_correct_yz(struct cuda_mfields* cmflds,
                            struct cuda_mfields* cmf, int p, Float3 fac,
                            Int3 ly, Int3 ry, Int3 lz, Int3 rz)
{
#if 0
  cuda_marder_correct_yz_gold(mflds, mf, p, fac, ly, ry, lz, rz);
  return;
#endif

  if (cmflds->n_patches() == 0) {
    return;
  }

  int my = cmflds->im(1);
  int mz = cmflds->im(2);

  int grid[2] = {(my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1]);

  marder_correct_yz<<<dimGrid, dimBlock>>>(
    cmflds->storage().view(_all, _all, _all, _all, p).to_kernel(),
    cmf->storage().view(_all, _all, _all, 0, p).to_kernel(), fac[1], fac[2],
    ly[1], ly[2], ry[1], ry[2], lz[1], lz[2], rz[1], rz[2], my, mz);
  cuda_sync_if_enabled();
}
