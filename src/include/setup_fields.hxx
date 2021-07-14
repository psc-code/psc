
#pragma once

#include "fields.hxx"
#include "centering.hxx"
#include <algorithm>

// ======================================================================
// SetupFields

namespace detail
{

template <typename MF>
struct SetupFields
{
  using Mfields = MF;
  using real_t = typename Mfields::real_t;

  template <typename FUNC>
  static void run(Mfields& mf, FUNC&& func)
  {
    const auto& grid = mf.grid();
    mpi_printf(grid.comm(), "**** Setting up fields...\n");

    for (int p = 0; p < mf.n_patches(); ++p) {
      auto& patch = grid.patches[p];
      auto F = make_Fields3d<dim_xyz>(mf[p]);

      int n_ghosts =
        std::max({mf.ibn()[0], mf.ibn()[1], mf.ibn()[2]}); // FIXME, not pretty
      // FIXME, do we need the ghost points?
      grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
        Int3 index{jx, jy, jz};

        for (int c = 0; c < 3; c++) {
          F(HX + c, jx, jy, jz) +=
            func(HX + c, Centering::getPos(patch, index, Centering::FC, c));
          F(EX + c, jx, jy, jz) +=
            func(EX + c, Centering::getPos(patch, index, Centering::EC, c));
          F(JXI + c, jx, jy, jz) +=
            func(JXI + c, Centering::getPos(patch, index, Centering::EC, c));
        }
      });
    }
  }

  template <typename FUNC>
  static void runScalar(Mfields& mf, FUNC&& func,
                        const Centering::Centerer& centerer)
  {
    const auto& grid = mf.grid();
    mpi_printf(grid.comm(), "**** Setting up scalar field...\n");

    for (int p = 0; p < mf.n_patches(); ++p) {
      auto& patch = grid.patches[p];
      auto F = make_Fields3d<dim_xyz>(mf[p]);

      int n_ghosts =
        std::max({mf.ibn()[0], mf.ibn()[1], mf.ibn()[2]}); // FIXME, not pretty
      // FIXME, do we need the ghost points?
      grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
        F(0, jx, jy, jz) += func(0, centerer.getPos(patch, {jx, jy, jz}));
      });
    }
  }
};

} // namespace detail

template <typename MF, typename FUNC>
void setupFields(MF& mflds, FUNC&& func)
{
  detail::SetupFields<MF>::run(mflds, std::forward<FUNC>(func));
}

template <typename MF, typename FUNC>
void setupScalarField(MF& fld, const Centering::Centerer& centerer, FUNC&& func)
{
  detail::SetupFields<MF>::runScalar(fld, std::forward<FUNC>(func), centerer);
}

#ifdef USE_CUDA
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#endif
