
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "fields_item_dive_cuda.hxx"

#define BND 2
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

// ======================================================================

void cuda_mfields_calc_dive_yz(MfieldsStateCuda& mflds, MfieldsCuda& mdive)
{
  auto dx = mflds.grid().domain.dx;

  assert(mflds.ibn() == Int3({0, BND, BND}));
  for (int p = 0; p < mflds.grid().n_patches(); p++) {
    auto flds = mflds.gt().view(_all, _all, _all, _all, p).to_kernel();
    auto bd = mdive.ibn();
    auto dive =
      mdive.gt()
        .view(_s(bd[0], -bd[0]), _s(bd[1], -bd[1]), _s(bd[2], -bd[2]), 0, p)
        .to_kernel();
    gt::launch<3>(dive.shape(), GT_LAMBDA(int i, int j, int k) {
      int iy = j + BND;
      int iz = k + BND;

      dive(0, j, k) = ((flds(0, iy, iz, EY) - flds(0, iy - 1, iz, EY)) / dx[1] +
                       (flds(0, iy, iz, EZ) - flds(0, iy, iz - 1, EZ)) / dx[2]);
    });
  }
  cuda_sync_if_enabled();
}
