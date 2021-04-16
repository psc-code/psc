
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

void cuda_marder_correct_xyz(MfieldsStateCuda& mflds, MfieldsCuda& mf, int p,
                             Float3 fac, Int3 lx, Int3 rx, Int3 ly, Int3 ry,
                             Int3 lz, Int3 rz)
{
  if (mflds.n_patches() == 0) {
    return;
  }

  auto gt_flds = mflds.gt().view(_all, _all, _all, _all, p).to_kernel();
  auto gt_f = mf.gt().view(_all, _all, _all, 0, p).to_kernel();

  gt::launch<3>(
    {gt_flds.shape(0), gt_flds.shape(1), gt_flds.shape(2)},
    GT_LAMBDA(int ix, int iy, int iz) {
      if ((ix - BND >= -lx[0] && ix - BND < rx[0]) &&
          (iy - BND >= -lx[1] && iy - BND < rx[1]) &&
          (iz - BND >= -lx[2] && iz - BND < rx[2])) {
        gt_flds(ix, iy, iz, EX) =
          gt_flds(ix, iy, iz, EX) +
          fac[0] * (gt_f(ix, iy + 1, iz) - gt_f(ix, iy, iz));
      }

      if ((ix - BND >= -ly[0] && ix - BND < ry[0]) &&
          (iy - BND >= -ly[1] && iy - BND < ry[1]) &&
          (iz - BND >= -ly[2] && iz - BND < ry[2])) {
        gt_flds(ix, iy, iz, EY) =
          gt_flds(ix, iy, iz, EY) +
          fac[1] * (gt_f(ix, iy + 1, iz) - gt_f(ix, iy, iz));
      }

      if ((ix - BND >= -lz[0] && ix - BND < rz[0]) &&
          (iy - BND >= -lz[1] && iy - BND < rz[1]) &&
          (iz - BND >= -lz[2] && iz - BND < rz[2])) {
        gt_flds(ix, iy, iz, EZ) =
          gt_flds(ix, iy, iz, EZ) +
          fac[2] * (gt_f(ix, iy, iz + 1) - gt_f(ix, iy, iz));
      }
    });
  cuda_sync_if_enabled();
}
