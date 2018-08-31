
#include "vpic_iface.h"

#include <psc_particles_single.h>
#include <psc_particles_single_by_kind.h>


// ======================================================================
// conversion

template<typename MP>
struct ConvertVpic;

template<typename MP>
struct ConvertToVpic;

template<typename MP>
struct ConvertFromVpic;

template<typename MP>
static void copy_from(MparticlesVpic& mprts, MP& mprts_from);

// ----------------------------------------------------------------------
// conversion to "single"

template<>
struct ConvertVpic<MparticlesSingle>
{
  ConvertVpic(MparticlesSingle& mprts_other, Grid& grid, int p)
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
  
  void operator()(struct vpic_mparticles_prt *prt, int n)
  {
    auto& prts_other = mprts_other_[p_];
    auto& prt_other = prts_other[n];
    
    assert(prt_other.kind < mprts_other_.grid().kinds.size());
    int i3[3];
    for (int d = 0; d < 3; d++) {
      float val = prt_other.x[d] / dx[d];
      i3[d] = (int) val;
      //mprintf("d %d val %g xi %g\n", d, val, prt_other.xi);
      assert(i3[d] >= -1 && i3[d] < im[d] + 1);
      prt->dx[d] = (val - i3[d]) * 2.f - 1.f;
      i3[d] += 1;
    }
    prt->i     = (i3[2] * im[1] + i3[1]) * im[0] + i3[0];
    prt->ux[0] = prt_other.p[0];
    prt->ux[1] = prt_other.p[1];
    prt->ux[2] = prt_other.p[2];
    prt->w     = prts_other.prt_wni(prt_other) / dVi;
    prt->kind  = prt_other.kind;
  }
};

template<>
struct ConvertFromVpic<MparticlesSingle> : ConvertVpic<MparticlesSingle>
{
  using Base = ConvertVpic<MparticlesSingle>;

  using Base::Base;
  
  void operator()(struct vpic_mparticles_prt *prt, int n)
  {
    auto& grid = mprts_other_.grid();
    auto& prts_other = mprts_other_[p_];
    
    assert(prt->kind < mprts_other_.grid().kinds.size());
    int i = prt->i;
    int i3[3];
    i3[2] = i / (im[0] * im[1]); i -= i3[2] * (im[0] * im[1]);
    i3[1] = i / im[0]; i-= i3[1] * im[0];
    i3[0] = i;
    auto xi = Vec3<float>{(i3[0] - 1 + .5f * (1.f + prt->dx[0])) * dx[0],
			  (i3[1] - 1 + .5f * (1.f + prt->dx[1])) * dx[1],
			  (i3[2] - 1 + .5f * (1.f + prt->dx[2])) * dx[2]};
    auto pxi = Vec3<float>{prt->ux[0], prt->ux[1], prt->ux[2]};
    auto kind = prt->kind;
    auto wni = float(prt->w * dVi);
    prts_other[n] = MparticlesSingle::particle_t{xi, pxi, wni, kind};
  }
};

#if 0
// ----------------------------------------------------------------------
// conversion to "single_by_kind"

template<>
struct ConvertVpic<MparticlesSingleByKind>
{
  ConvertVpic(MparticlesSingleByKind& mprts_other, Grid& grid, int p)
    : mprts_other_(mprts_other), p_(p)
  {
  }

protected:
  MparticlesSingleByKind& mprts_other_;
  int p_;
};

template<>
struct ConvertToVpic<MparticlesSingleByKind> : ConvertVpic<MparticlesSingleByKind>
{
  using Base = ConvertVpic<MparticlesSingleByKind>;

  using Base::Base;
  
  void operator()(struct vpic_mparticles_prt *prt, int n)
  {
    particle_single_by_kind *part = &mprts_other_.bkmprts->at(p_, n);
    
    //  assert(part->kind < ppsc->nr_kinds);
    prt->dx[0] = part->dx[0];
    prt->dx[1] = part->dx[1];
    prt->dx[2] = part->dx[2];
    prt->i     = part->i;
    prt->ux[0] = part->ux[0];
    prt->ux[1] = part->ux[1];
    prt->ux[2] = part->ux[2];
    prt->w     = part->w;
    prt->kind  = part->kind;
  }
};

template<>
struct ConvertFromVpic<MparticlesSingleByKind> : ConvertVpic<MparticlesSingleByKind>
{
  using Base = ConvertVpic<MparticlesSingleByKind>;

  using Base::Base;
  
  void operator()(struct vpic_mparticles_prt *prt, int n)
  {
    particle_single_by_kind *part = &mprts_other_.bkmprts->at(p_, n);
    
    assert(prt->kind < 2); // FIXMEppsc->nr_kinds);
    part->dx[0] = prt->dx[0];
    part->dx[1] = prt->dx[1];
    part->dx[2] = prt->dx[2];
    part->i     = prt->i;
    part->ux[0] = prt->ux[0];
    part->ux[1] = prt->ux[1];
    part->ux[2] = prt->ux[2];
    part->w     = prt->w;
    part->kind  = prt->kind;
  }
};
#endif
  
template<typename MP>
void copy_to(MparticlesVpic& mprts_from, MP& mprts_to)
{
  Particles& vmprts = mprts_from.vmprts();
  int n_patches = mprts_to.n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from.get_size_all(n_prts_by_patch);
  mprts_to.reserve_all(n_prts_by_patch);
  mprts_to.resize_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    ConvertFromVpic<MP> convert_from_vpic(mprts_to, *vmprts.grid(), p);
    vpic_mparticles_get_particles(&vmprts, n_prts, off, convert_from_vpic);

    off += n_prts;
  }
}

template<>
void copy_from(MparticlesVpic& mprts_to, MparticlesSingle& mprts_from)
{
  Particles& vmprts = mprts_to.vmprts();
  int n_patches = mprts_to.n_patches();
  uint n_prts_by_patch[n_patches];
  // reset particle counts to zero, then use push_back to add back new particles
  for (int p = 0; p < n_patches; p++) {
    n_prts_by_patch[p] = 0;
  }
  mprts_to.resize_all(n_prts_by_patch);

  mprts_from.get_size_all(n_prts_by_patch);
  mprts_to.reserve_all(n_prts_by_patch);
  
  for (int p = 0; p < n_patches; p++) {
    ConvertToVpic<MparticlesSingle> convert_to_vpic(mprts_from, *vmprts.grid(), p);

    int n_prts = n_prts_by_patch[p];
    for (int n = 0; n < n_prts; n++) {
      struct vpic_mparticles_prt prt;
      convert_to_vpic(&prt, n);
      mprts_to.push_back(&prt);
    }
  }
}

#if 0
template<>
void copy_from(MparticlesVpic& mprts_to, MparticlesSingleByKind& mprts_from)
{
  Particles& vmprts = mprts_to.vmprts();
  int n_patches = mprts_to.n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from.get_size_all(n_prts_by_patch);
  mprts_to.reserve_all(n_prts_by_patch);
  mprts_to.resize_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    ConvertToVpic<MparticlesSingleByKind> convert_to_vpic(mprts_from, *vmprts.grid(), p);
    vpic_mparticles_set_particles(&vmprts, n_prts, off, convert_to_vpic);

    off += n_prts;
  }
}
#endif

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
  //  { std::type_index(typeid(MparticlesSingleByKind)), psc_mparticles_vpic_copy_to<MparticlesSingleByKind> },
};

template<>
const MparticlesVpic::Convert MparticlesVpic::convert_from_ = {
  { std::type_index(typeid(MparticlesSingle)), psc_mparticles_vpic_copy_from<MparticlesSingle> },
  //  { std::type_index(typeid(MparticlesSingleByKind)), psc_mparticles_vpic_copy_from<MparticlesSingleByKind> },
};

