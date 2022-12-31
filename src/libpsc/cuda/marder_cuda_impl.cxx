
#include "marder_cuda_impl.hxx"

#define BND 2

void cuda_marder_correct_yz(MfieldsStateCuda::Storage& gt_mflds,
                            MfieldsCuda::Storage& gt_mf, int p, Float3 fac,
                            Int3 ly, Int3 ry, Int3 lz, Int3 rz)
{
  auto gt_flds = gt_mflds.view(_all, _all, _all, _all, p).to_kernel();
  auto gt_f = gt_mf.view(_all, _all, _all, 0, p).to_kernel();
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

void cuda_marder_correct_xyz(MfieldsStateCuda::Storage& gt_mflds,
                             MfieldsCuda::Storage& gt_mf, int p, Float3 fac,
                             Int3 lx, Int3 rx, Int3 ly, Int3 ry, Int3 lz,
                             Int3 rz)
{
  auto gt_flds = gt_mflds.view(_all, _all, _all, _all, p).to_kernel();
  auto gt_f = gt_mf.view(_all, _all, _all, 0, p).to_kernel();
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
