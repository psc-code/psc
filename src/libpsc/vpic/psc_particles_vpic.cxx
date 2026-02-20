
#include "vpic_iface.h"

#include <psc_particles_single.h>

#ifdef DO_VPIC
using VpicConfig = VpicConfigWrap;
#else
using VpicConfig = VpicConfigPsc;
#endif

using Grid = VpicConfig::Grid;
using MparticlesVpic = VpicConfig::Mparticles;

// ======================================================================
// conversion

template <typename Mparticles>
void copy_to(MparticlesBase& mprts_from_base, MparticlesBase& mprts_to_base)
{
  auto& mprts_from = dynamic_cast<MparticlesVpic&>(mprts_from_base);
  auto& mprts_to = dynamic_cast<Mparticles&>(mprts_to_base);
  auto& vgrid = mprts_from.vgrid();
  int im[3] = {vgrid.nx + 2, vgrid.ny + 2, vgrid.nz + 2};
  float dx[3] = {vgrid.dx, vgrid.dy, vgrid.dz};
  float dVi = 1.f / (dx[0] * dx[1] * dx[2]); // FIXME, vgrid->dVi?

  auto n_prts_by_patch = mprts_from.sizeByPatch();
  mprts_to.reserve_all(n_prts_by_patch);
  mprts_to.clear();

  for (int p = 0; p < mprts_to.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];

    auto prts_from = mprts_from[p];
    for (auto sp = prts_from.cbegin(); sp != prts_from.cend(); ++sp) {
      assert(sp->id < mprts_to.grid().kinds.size());

      for (unsigned int n = 0; n < sp->np; n++) {
        auto& vprt = sp->p[n];
        int i = vprt.i;
        int i3[3];
        i3[2] = i / (im[0] * im[1]);
        i -= i3[2] * (im[0] * im[1]);
        i3[1] = i / im[0];
        i -= i3[1] * im[0];
        i3[0] = i;
        auto x = Vec3<float>{(i3[0] - 1 + .5f * (1.f + vprt.dx)) * dx[0],
                             (i3[1] - 1 + .5f * (1.f + vprt.dy)) * dx[1],
                             (i3[2] - 1 + .5f * (1.f + vprt.dz)) * dx[2]};
        auto u = Vec3<float>{vprt.ux, vprt.uy, vprt.uz};
        auto kind = sp->id;
        auto qni_wni =
          float(vprt.w * dVi) * float(mprts_to.grid().kinds[kind].q);
        mprts_to.push_back(
          p, {x, u, qni_wni, kind, psc::particle::Id{}, psc::particle::Tag{}});
      }
    }
  }
}

template <typename Mparticles>
void copy_from(MparticlesBase& mprts_to_base, MparticlesBase& mprts_from_base)
{
  auto& mprts_to = dynamic_cast<MparticlesVpic&>(mprts_to_base);
  auto& mprts_from = dynamic_cast<Mparticles&>(mprts_from_base);
  auto& vgrid = mprts_to.vgrid();
  int im[3] = {vgrid.nx + 2, vgrid.ny + 2, vgrid.nz + 2};
  float dx[3] = {vgrid.dx, vgrid.dy, vgrid.dz};
  float dVi = 1.f / (dx[0] * dx[1] * dx[2]); // FIXME, vgrid->dVi?

  mprts_to.reset();
  mprts_to.reserve_all(mprts_from.sizeByPatch());

  auto accessor_from = mprts_from.accessor();
  for (int p = 0; p < mprts_to.n_patches(); p++) {
    auto prts_from = accessor_from[p];

    int n_prts = prts_from.size();
    for (auto prt_from : prts_from) {
      assert(prt_from.kind() < mprts_to.grid().kinds.size());

      int i3[3];
      float dx[3];
      for (int d = 0; d < 3; d++) {
        float val = prt_from.x()[d] / dx[d];
        i3[d] = (int)val;
        // mprintf("d %d val %g xi %g\n", d, val, prt_from.xi);
        assert(i3[d] >= -1 && i3[d] < im[d] + 1);
        dx[d] = (val - i3[d]) * 2.f - 1.f;
        i3[d] += 1;
      }
      typename MparticlesVpic::Particle prt;
      prt.dx = dx[0];
      prt.dy = dx[1];
      prt.dz = dx[2];
      prt.i = (i3[2] * im[1] + i3[1]) * im[0] + i3[0];
      prt.ux = prt_from.u()[0];
      prt.uy = prt_from.u()[1];
      prt.uz = prt_from.u()[2];
      prt.w = prt_from.w() / dVi;
      mprts_to.push_back(prt_from.kind(), prt);
    }
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_methods

template <>
const MparticlesVpic::Convert MparticlesVpic::convert_to_ = {
  {std::type_index(typeid(MparticlesSingle)), copy_to<MparticlesSingle>},
};

template <>
const MparticlesVpic::Convert MparticlesVpic::convert_from_ = {
  {std::type_index(typeid(MparticlesSingle)), copy_from<MparticlesSingle>},
};
