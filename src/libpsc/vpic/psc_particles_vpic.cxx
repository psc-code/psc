
#include "psc_particles_vpic.h"

#include <psc_particles_single_by_kind.h>
#include <psc_method.h>

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// vpic_mparticles_get_size_all

void vpic_mparticles_get_size_all(Particles *vmprts, int n_patches,
				  uint *n_prts_by_patch)
{
  assert(n_patches == 1);
  uint n_prts = 0;

  for (auto sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    n_prts += sp->np;
  }

  n_prts_by_patch[0] = n_prts;
}

// ======================================================================
// conversion

template<typename MP>
struct ConvertVpic;

template<typename MP>
struct ConvertToVpic;

template<typename MP>
struct ConvertFromVpic;

template<typename MP>
static void copy_from(mparticles_vpic_t mprts, MP mprts_from);

// ----------------------------------------------------------------------
// conversion to "single"

template<>
struct ConvertVpic<mparticles_single_t>
{
  ConvertVpic(mparticles_single_t& mprts_other, Grid& grid, int p)
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
  mparticles_single_t& mprts_other_;
  int p_;
  int im[3];
  float dx[3];
  float dVi;
};

template<>
struct ConvertToVpic<mparticles_single_t> : ConvertVpic<mparticles_single_t>
{
  using Base = ConvertVpic<mparticles_single_t>;

  using Base::Base;
  
  void operator()(struct vpic_mparticles_prt *prt, int n)
  {
    particle_single_t *part = &mprts_other_[p_][n];
    
    assert(part->kind() < ppsc->nr_kinds);
    int i3[3];
    for (int d = 0; d < 3; d++) {
      float val = (&part->xi)[d] / dx[d];
      i3[d] = (int) val;
      //mprintf("d %d val %g xi %g\n", d, val, part->xi);
      assert(i3[d] >= -1 && i3[d] < im[d] + 1);
      prt->dx[d] = (val - i3[d]) * 2.f - 1.f;
      i3[d] += 1;
    }
    prt->i     = (i3[2] * im[1] + i3[1]) * im[0] + i3[0];
    prt->ux[0] = part->pxi;
    prt->ux[1] = part->pyi;
    prt->ux[2] = part->pzi;
    prt->w     = part->qni_wni_ / ppsc->kinds[part->kind_].q / dVi;
    prt->kind  = part->kind();
  }
};

template<>
struct ConvertFromVpic<mparticles_single_t> : ConvertVpic<mparticles_single_t>
{
  using Base = ConvertVpic<mparticles_single_t>;

  using Base::Base;
  
  void operator()(struct vpic_mparticles_prt *prt, int n)
  {
    particle_single_t *part = &mprts_other_[p_][n];
    
    assert(prt->kind < ppsc->nr_kinds);
    int i = prt->i;
    int i3[3];
    i3[2] = i / (im[0] * im[1]); i -= i3[2] * (im[0] * im[1]);
    i3[1] = i / im[0]; i-= i3[1] * im[0];
    i3[0] = i;
    part->xi      = (i3[0] - 1 + .5f * (1.f + prt->dx[0])) * dx[0];
    part->yi      = (i3[1] - 1 + .5f * (1.f + prt->dx[1])) * dx[1];
    part->zi      = (i3[2] - 1 + .5f * (1.f + prt->dx[2])) * dx[2];
    float w = part->zi / dx[2];
    if (!(w >= 0 && w <= im[2] - 2)) {
      printf("w %g im %d i3 %d dx %g\n", w, im[2], i3[2], prt->dx[2]);
    }
    assert(w >= 0 && w <= im[2] - 2);
    part->kind_   = prt->kind;
    part->pxi     = prt->ux[0];
    part->pyi     = prt->ux[1];
    part->pzi     = prt->ux[2];
    part->qni_wni_ = ppsc->kinds[prt->kind].q * prt->w * dVi;
  }
};

// ----------------------------------------------------------------------
// conversion to "single_by_kind"

template<>
struct ConvertVpic<mparticles_single_by_kind_t>
{
  ConvertVpic(mparticles_single_by_kind_t& mprts_other, Grid& grid, int p)
    : mprts_other_(mprts_other), p_(p)
  {
  }

protected:
  mparticles_single_by_kind_t& mprts_other_;
  int p_;
};

template<>
struct ConvertToVpic<mparticles_single_by_kind_t> : ConvertVpic<mparticles_single_by_kind_t>
{
  using Base = ConvertVpic<mparticles_single_by_kind_t>;

  using Base::Base;
  
  void operator()(struct vpic_mparticles_prt *prt, int n)
  {
    particle_single_by_kind *part = &mprts_other_->bkmprts->at(p_, n);
    
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
struct ConvertFromVpic<mparticles_single_by_kind_t> : ConvertVpic<mparticles_single_by_kind_t>
{
  using Base = ConvertVpic<mparticles_single_by_kind_t>;

  using Base::Base;
  
  void operator()(struct vpic_mparticles_prt *prt, int n)
  {
    particle_single_by_kind *part = &mprts_other_->bkmprts->at(p_, n);
    
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
  
template<typename MP>
void copy_to(mparticles_vpic_t mprts_from, MP mprts_to)
{
  Particles *vmprts = mprts_from->vmprts;
  int n_patches = mprts_to->n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from->get_size_all(n_prts_by_patch);
  mprts_to->reserve_all(n_prts_by_patch);
  mprts_to->resize_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    ConvertFromVpic<MP> convert_from_vpic(mprts_to, *vmprts->grid(), p);
    vpic_mparticles_get_particles(vmprts, n_prts, off, convert_from_vpic);

    off += n_prts;
  }
}

template<>
void copy_from<mparticles_single_t>(mparticles_vpic_t mprts_to, mparticles_single_t mprts_from)
{
  Particles *vmprts = mprts_to->vmprts;
  int n_patches = mprts_to->n_patches();
  uint n_prts_by_patch[n_patches];
  // reset particle counts to zero, then use push_back to add back new particles
  for (int p = 0; p < n_patches; p++) {
    n_prts_by_patch[p] = 0;
  }
  mprts_to->resize_all(n_prts_by_patch);

  mprts_from->get_size_all(n_prts_by_patch);
  mprts_to->reserve_all(n_prts_by_patch);
  
  for (int p = 0; p < n_patches; p++) {
    ConvertToVpic<mparticles_single_t> convert_to_vpic(mprts_from, *vmprts->grid(), p);

    int n_prts = n_prts_by_patch[p];
    for (int n = 0; n < n_prts; n++) {
      struct vpic_mparticles_prt prt;
      convert_to_vpic(&prt, n);
      Simulation_mprts_push_back(mprts_to->sim, vmprts, &prt);
    }
  }
}

template<>
void copy_from<mparticles_single_by_kind_t>(mparticles_vpic_t mprts_to, mparticles_single_by_kind_t mprts_from)
{
  Particles *vmprts = mprts_to->vmprts;
  int n_patches = mprts_to->n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from->get_size_all(n_prts_by_patch);
  mprts_to->reserve_all(n_prts_by_patch);
  mprts_to->resize_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    ConvertToVpic<mparticles_single_by_kind_t> convert_to_vpic(mprts_from, *vmprts->grid(), p);
    vpic_mparticles_set_particles(vmprts, n_prts, off, convert_to_vpic);

    off += n_prts;
  }
}

template<typename MP>
static void psc_mparticles_vpic_copy_from(struct psc_mparticles *mprts,
					  struct psc_mparticles *mprts_other,
					  unsigned int flags)
{
  copy_from(mparticles_vpic_t(mprts), MP(mprts_other));
}

template<typename MP>
static void psc_mparticles_vpic_copy_to(struct psc_mparticles *mprts,
					struct psc_mparticles *mprts_other,
					unsigned int flags)
{
  copy_to(mparticles_vpic_t(mprts), MP(mprts_other));
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_methods

static struct mrc_obj_method psc_mparticles_vpic_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"          , psc_mparticles_vpic_copy_to<mparticles_single_t>),
  MRC_OBJ_METHOD("copy_from_single"        , psc_mparticles_vpic_copy_from<mparticles_single_t>),
  MRC_OBJ_METHOD("copy_to_single_by_kind"  , psc_mparticles_vpic_copy_to<mparticles_single_by_kind_t>),
  MRC_OBJ_METHOD("copy_from_single_by_kind", psc_mparticles_vpic_copy_from<mparticles_single_by_kind_t>),
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_vpic_setup

// FIXME
#define DONTEXPAND

static void
psc_mparticles_vpic_setup(struct psc_mparticles *_mprts)
{
  mparticles_vpic_t mprts(_mprts);

  new(mprts.sub()) psc_mparticles_vpic{ppsc->grid};

  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &mprts->sim);
  mprts->vmprts = Simulation_get_particles(mprts->sim);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_get_size_all

static void
psc_mparticles_vpic_get_size_all(struct psc_mparticles *mprts, uint *n_prts_by_patch)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);
  
  vpic_mparticles_get_size_all(sub->vmprts, mprts->nr_patches, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_reserve_all

static void
psc_mparticles_vpic_reserve_all(struct psc_mparticles *mprts, uint *n_prts_by_patch)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);

  Simulation_mprts_reserve_all(sub->sim, sub->vmprts, mprts->nr_patches, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_resize_all

static void
psc_mparticles_vpic_resize_all(struct psc_mparticles *mprts, uint *n_prts_by_patch)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);

  Simulation_mprts_resize_all(sub->sim, sub->vmprts, mprts->nr_patches, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_get_nr_particles

static unsigned int
psc_mparticles_vpic_get_nr_particles(struct psc_mparticles *mprts)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);
  
  return Simulation_mprts_get_nr_particles(sub->sim, sub->vmprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_inject

static void
psc_mparticles_vpic_inject(struct psc_mparticles *mprts, int p,
			   const struct psc_particle_inject *prt)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);

  Simulation_inject_particle(sub->sim, sub->vmprts, p, prt);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "vpic"
  
struct psc_mparticles_ops_vpic : psc_mparticles_ops {
  psc_mparticles_ops_vpic() {
    name                    = "vpic";
    size                    = sizeof(struct psc_mparticles_vpic);
    methods                 = psc_mparticles_vpic_methods;
    setup                   = psc_mparticles_vpic_setup;
    reserve_all             = psc_mparticles_vpic_reserve_all;
    get_size_all            = psc_mparticles_vpic_get_size_all;
    resize_all              = psc_mparticles_vpic_resize_all;
    get_nr_particles        = psc_mparticles_vpic_get_nr_particles;
    inject                  = psc_mparticles_vpic_inject;
  }
} psc_mparticles_vpic_ops;




