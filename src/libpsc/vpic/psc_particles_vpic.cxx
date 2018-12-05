
#include "vpic_iface.h"

#include <psc_particles_single.h>
#include <psc_particles_single_by_kind.h>

#ifdef DO_VPIC
using VpicConfig = VpicConfigWrap;
#else
using VpicConfig = VpicConfigPsc;
#endif

using Grid = VpicConfig::Grid;
using MparticlesVpic = VpicConfig::Mparticles;

// ======================================================================
// conversion

template<typename MP>
struct ConvertVpic;

template<typename MP>
struct ConvertToVpic;

template<typename MP>
static void copy_from(MparticlesVpic& mprts, MP& mprts_from);

// ----------------------------------------------------------------------
// conversion to "single"

template<>
struct ConvertVpic<MparticlesSingle>
{
  ConvertVpic(MparticlesSingle& mprts_other, const Grid& grid, int p)
    : mprts_other_(mprts_other), p_(p)
  {
    im[0] = grid.nx + 2;
    im[1] = grid.ny + 2;
    im[2] = grid.nz + 2;
    dx[0] = grid.dx;
    dx[1] = grid.dy;
    dx[2] = grid.dz;
    dVi = 1.f / (dx[0] * dx[1] * dx[2]);
  }

protected:
  MparticlesSingle& mprts_other_;
  int p_;
  int im[3];
  float dx[3];
  float dVi;
};

template<>
struct ConvertToVpic<MparticlesSingle> : ConvertVpic<MparticlesSingle>
{
  using Base = ConvertVpic<MparticlesSingle>;

  using Base::Base;
  
  void operator()(vpic_mparticles_prt& prt, int n)
  {
    auto& prts_other = mprts_other_[p_];
    auto& prt_other = prts_other[n];
    
    assert(prt_other.kind() < mprts_other_.grid().kinds.size());
    int i3[3];
    for (int d = 0; d < 3; d++) {
      float val = prt_other.x()[d] / dx[d];
      i3[d] = (int) val;
      //mprintf("d %d val %g xi %g\n", d, val, prt_other.xi);
      assert(i3[d] >= -1 && i3[d] < im[d] + 1);
      prt.dx[d] = (val - i3[d]) * 2.f - 1.f;
      i3[d] += 1;
    }
    prt.i     = (i3[2] * im[1] + i3[1]) * im[0] + i3[0];
    prt.ux[0] = prt_other.u()[0];
    prt.ux[1] = prt_other.u()[1];
    prt.ux[2] = prt_other.u()[2];
    prt.w     = prts_other.prt_wni(prt_other) / dVi;
    prt.kind  = prt_other.kind();
  }
};

template<typename MP>
void copy_to(MparticlesVpic& mprts_from, MP& mprts_to)
{
  auto& vgrid = mprts_from.vgrid();
  int im[3] = { vgrid.nx + 2, vgrid.ny + 2, vgrid.nz + 2 };
  float dx[3] = { vgrid.dx, vgrid.dy, vgrid.dz };
  float dVi = 1.f / (dx[0] * dx[1] * dx[2]); // FIXME, vgrid->dVi?

  auto n_prts_by_patch = mprts_from.get_size_all();
  mprts_to.reserve_all(n_prts_by_patch);
  mprts_to.reset();
  
  for (int p = 0; p < mprts_to.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];

    for (auto sp = mprts_from.cbegin(); sp != mprts_from.cend(); ++sp) {
      assert(sp->id < mprts_to.grid().kinds.size());

      for (unsigned int n = 0; n < sp->np; n++) {
	auto& vprt = sp->p[n];
	int i = vprt.i;
	int i3[3];
	i3[2] = i / (im[0] * im[1]); i -= i3[2] * (im[0] * im[1]);
	i3[1] = i / im[0]; i-= i3[1] * im[0];
	i3[0] = i;
	auto x = Vec3<float>{(i3[0] - 1 + .5f * (1.f + vprt.dx)) * dx[0],
			     (i3[1] - 1 + .5f * (1.f + vprt.dy)) * dx[1],
			     (i3[2] - 1 + .5f * (1.f + vprt.dz)) * dx[2]};
	auto u = Vec3<float>{vprt.ux, vprt.uy, vprt.uz};
	auto kind = sp->id;
	auto qni_wni = float(vprt.w * dVi) * float(mprts_to.grid().kinds[kind].q);
	mprts_to[p].push_back({x, u, qni_wni, kind});
      }
    }
  }
}

template<>
void copy_from(MparticlesVpic& mprts_to, MparticlesSingle& mprts_from)
{
  auto& vgrid = mprts_to.vgrid();
  int im[3] = { vgrid.nx + 2, vgrid.ny + 2, vgrid.nz + 2 };
  float dx[3] = { vgrid.dx, vgrid.dy, vgrid.dz };
  float dVi = 1.f / (dx[0] * dx[1] * dx[2]); // FIXME, vgrid->dVi?

  mprts_to.reset();

  auto n_prts_by_patch = mprts_from.get_size_all();
  mprts_to.reserve_all(n_prts_by_patch);
  
  for (int p = 0; p < mprts_to.n_patches(); p++) {
    auto& prts_from = mprts_from[p];

    int n_prts = n_prts_by_patch[p];
    for (int n = 0; n < n_prts; n++) {
      auto& prt_from = prts_from[n];
      assert(prt_from.kind() < mprts_to.grid().kinds.size());

      vpic_mparticles_prt prt;
      int i3[3];
      for (int d = 0; d < 3; d++) {
	float val = prt_from.x()[d] / dx[d];
	i3[d] = (int) val;
	//mprintf("d %d val %g xi %g\n", d, val, prt_from.xi);
	assert(i3[d] >= -1 && i3[d] < im[d] + 1);
	prt.dx[d] = (val - i3[d]) * 2.f - 1.f;
	i3[d] += 1;
      }
      prt.i     = (i3[2] * im[1] + i3[1]) * im[0] + i3[0];
      prt.ux[0] = prt_from.u()[0];
      prt.ux[1] = prt_from.u()[1];
      prt.ux[2] = prt_from.u()[2];
      prt.w     = prts_from.prt_wni(prt_from) / dVi;
      prt.kind  = prt_from.kind();
      mprts_to.push_back(prt);
    }
  }
}

template<typename MP>
static void psc_mparticles_vpic_copy_from(MparticlesBase& mp,
					  MparticlesBase& mp_other)
{
  copy_from(dynamic_cast<MparticlesVpic&>(mp), dynamic_cast<MP&>(mp_other));
}

template<typename MP>
static void psc_mparticles_vpic_copy_to(MparticlesBase& mp,
					MparticlesBase& mp_other)
{
  copy_to(dynamic_cast<MparticlesVpic&>(mp), dynamic_cast<MP&>(mp_other));
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_methods

template<>
const MparticlesVpic::Convert MparticlesVpic::convert_to_ = {
  { std::type_index(typeid(MparticlesSingle)), psc_mparticles_vpic_copy_to<MparticlesSingle> },
};

template<>
const MparticlesVpic::Convert MparticlesVpic::convert_from_ = {
  { std::type_index(typeid(MparticlesSingle)), psc_mparticles_vpic_copy_from<MparticlesSingle> },
};

