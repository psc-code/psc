
#include "psc_particles_vpic.h"

#include <psc_particles_single_by_kind.h>
#include <psc_method.h>

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// vpic_mparticles_get_size_all

void vpic_mparticles_get_size_all(Particles *vmprts, int n_patches,
				  int *n_prts_by_patch)
{
  assert(n_patches == 1);
  int n_prts = 0;

  for (auto sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    n_prts += sp->np;
  }

  n_prts_by_patch[0] = n_prts;
}

// ======================================================================
// conversion

struct copy_ctx
{
  copy_ctx(psc_mparticles *mprts, Grid& grid, int p)
    : mprts_(mprts), p_(p)
  {
    im[0] = grid.nx + 2;
    im[1] = grid.ny + 2;
    im[2] = grid.nz + 2;
    dx[0] = grid.dx;
    dx[1] = grid.dy;
    dx[2] = grid.dz;
    dVi = 1.f / (dx[0] * dx[1] * dx[2]);
  }
  
  struct psc_mparticles *mprts_;
  int p_;
  int im[3];
  float dx[3];
  float dVi;
};

struct copy2_ctx
{
  copy2_ctx(bk_mparticles *bkmprts, Grid& grid, int p)
    : bkmprts_(bkmprts), p_(p)
  {
  }
  
  bk_mparticles *bkmprts_;
  int p_;
};

template<typename MP, typename F>
static void copy_to(mparticles_vpic_t mprts, MP mprts_to, F convert_to)
{
  Particles *vmprts = mprts->vmprts;
  
  int n_prts_by_patch[mprts.n_patches()];
  vpic_mparticles_get_size_all(vmprts, mprts.n_patches(), n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx(mprts_to.mprts(), *vmprts->grid(), p);
    vpic_mparticles_get_particles(vmprts, n_prts, off, convert_to, &ctx);

    off += n_prts;
  }
}

template<typename MP, typename F>
static void copy2_to(mparticles_vpic_t mprts, MP mprts_to, F convert_to)
{
  Particles *vmprts = mprts->vmprts;
  bk_mparticles *bkmprts = mprts_to->bkmprts;

  assert(mprts.n_patches() == 1);
  int n_prts_by_patch[mprts.n_patches()];
  vpic_mparticles_get_size_all(vmprts, mprts.n_patches(), n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy2_ctx ctx(bkmprts, *vmprts->grid(), p);
    vpic_mparticles_get_particles(vmprts, n_prts, off, convert_to, &ctx);

    off += n_prts;
  }
}

template<typename MP, typename F>
static void copy_from(mparticles_vpic_t mprts, MP mprts_from, F convert_from)
{
  Particles *vmprts = mprts->vmprts;

  int n_prts_by_patch[mprts.n_patches()];
  for (int p = 0; p < mprts.n_patches(); p++) {
    n_prts_by_patch[p] = 0;
  }
  // reset particle counts to zero, then use push_back to add back new particles
  Simulation_mprts_resize_all(mprts->sim, vmprts, mprts.n_patches(), n_prts_by_patch);
  psc_mparticles_get_size_all(mprts_from.mprts(), n_prts_by_patch);
  
  for (int p = 0; p < mprts.n_patches(); p++) {
    struct copy_ctx ctx(mprts_from.mprts(), *vmprts->grid(), p);
    struct vpic_mparticles_prt prt;

    int n_prts = n_prts_by_patch[p];
    for (int n = 0; n < n_prts; n++) {
      convert_from(&prt, n, &ctx);
      Simulation_mprts_push_back(mprts->sim, vmprts, &prt);
    }
  }
}

template<typename MP, typename F>
static void copy2_from(mparticles_vpic_t mprts, MP mprts_from, F convert_from)
{
  Particles *vmprts = mprts->vmprts;
  bk_mparticles *bkmprts = mprts_from->bkmprts;

  int n_prts_by_patch[mprts.n_patches()];
  psc_mparticles_get_size_all(mprts_from.mprts(), n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy2_ctx ctx(bkmprts, *vmprts->grid(), p);
    vpic_mparticles_set_particles(vmprts, n_prts, off, convert_from, &ctx);

    off += n_prts;
  }
}

// ======================================================================
// conversion to "single"

struct ConvertFromSingle
{
  void operator()(struct vpic_mparticles_prt *prt, int n, void *_ctx)
  {
    struct copy_ctx *ctx = (struct copy_ctx *) _ctx;
    particle_single_t *part = &mparticles_single_t(ctx->mprts_)[ctx->p_][n];
    
    assert(part->kind() < ppsc->nr_kinds);
    int *im = ctx->im;
    float *dx = ctx->dx;
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
    prt->w     = part->qni_wni / ppsc->kinds[part->kind_].q / ctx->dVi;
    prt->kind  = part->kind();
  }
};

struct ConvertToSingle
{
  void operator()(struct vpic_mparticles_prt *prt, int n, void *_ctx)
  {
    struct copy_ctx *ctx = (struct copy_ctx *) _ctx;
    particle_single_t *part = &mparticles_single_t(ctx->mprts_)[ctx->p_][n];
    
    assert(prt->kind < ppsc->nr_kinds);
    int *im = ctx->im;
    float *dx = ctx->dx;
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
    part->qni_wni = ppsc->kinds[prt->kind].q * prt->w * ctx->dVi;
  }
};

static void
convert_from_single_by_kind(struct vpic_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy2_ctx *ctx = (struct copy2_ctx *) _ctx;
  particle_single_by_kind *part = &ctx->bkmprts_->at(ctx->p_, n);

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

static void
convert_to_single_by_kind(struct vpic_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy2_ctx *ctx = (struct copy2_ctx *) _ctx;
  particle_single_by_kind *part = &ctx->bkmprts_->at(ctx->p_, n);

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

static void
psc_mparticles_vpic_copy_from_single(struct psc_mparticles *mprts,
				     struct psc_mparticles *mprts_single, unsigned int flags)
{
  ConvertFromSingle convert_from_single;
  copy_from(mparticles_vpic_t(mprts), mparticles_single_t(mprts_single), convert_from_single);
}

static void
psc_mparticles_vpic_copy_to_single(struct psc_mparticles *mprts,
				   struct psc_mparticles *mprts_single, unsigned int flags)
{
  ConvertToSingle convert_to_single;
  copy_to(mparticles_vpic_t(mprts), mparticles_single_t(mprts_single), convert_to_single);
}

// ======================================================================
// conversion to "single_by_kind"

static void
psc_mparticles_vpic_copy_from_single_by_kind(struct psc_mparticles *mprts,
				    struct psc_mparticles *mprts_single_by_kind, unsigned int flags)
{
  copy2_from(mparticles_vpic_t(mprts), mparticles_single_by_kind_t(mprts_single_by_kind),
	     convert_from_single_by_kind);
}

static void
psc_mparticles_vpic_copy_to_single_by_kind(struct psc_mparticles *mprts,
				  struct psc_mparticles *mprts_single_by_kind, unsigned int flags)
{
  copy2_to(mparticles_vpic_t(mprts), mparticles_single_by_kind_t(mprts_single_by_kind),
	   convert_to_single_by_kind);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_methods

static struct mrc_obj_method psc_mparticles_vpic_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_mparticles_vpic_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mparticles_vpic_copy_from_single),
  MRC_OBJ_METHOD("copy_to_single_by_kind"  , psc_mparticles_vpic_copy_to_single_by_kind),
  MRC_OBJ_METHOD("copy_from_single_by_kind", psc_mparticles_vpic_copy_from_single_by_kind),
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_vpic_setup

static void
psc_mparticles_vpic_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);

  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sub->sim);
  sub->vmprts = Simulation_get_particles(sub->sim);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_get_size_all

static void
psc_mparticles_vpic_get_size_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);
  
  vpic_mparticles_get_size_all(sub->vmprts, mprts->nr_patches, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_reserve_all

static void
psc_mparticles_vpic_reserve_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);

  Simulation_mprts_reserve_all(sub->sim, sub->vmprts, mprts->nr_patches, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_resize_all

static void
psc_mparticles_vpic_resize_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
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




