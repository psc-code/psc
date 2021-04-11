
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

  auto gt_flds = cmflds->storage().view(_all, _all, _all, _all, p).to_kernel();
  auto gt_f = cmf->storage().view(_all, _all, _all, 0, p).to_kernel();

  gt::launch<2>({gt_flds.shape(1), gt_flds.shape(2)}, GT_LAMBDA(int iy,
                                                                int iz) {
    if (iy - BND >= -ly[1] && iy - BND < ry[1] && iz - BND >= -ly[2] &&
        iz - BND < ry[2]) {
      gt_flds(0, iy, iz, EY) = gt_flds(0, iy, iz, EY) +
                               fac[1] * (gt_f(0, iy + 1, iz) - gt_f(0, iy, iz));
    }

    if (iy - BND >= -lz[1] && iy - BND < rz[1] && iz - BND >= -lz[2] &&
        iz - BND < rz[2]) {
      gt_flds(0, iy, iz, EZ) = gt_flds(0, iy, iz, EZ) +
                               fac[2] * (gt_f(0, iy, iz + 1) - gt_f(0, iy, iz));
    }
  });
  cuda_sync_if_enabled();
}
