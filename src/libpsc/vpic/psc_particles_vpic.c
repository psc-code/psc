
#include "psc_particles_vpic.h"

#include "vpic_iface.h"

// ======================================================================
// conversion

struct copy_ctx {
  struct psc_mparticles *mprts;
  int p;
  int im[3];
  float dx[3];
  float dVi;
};

static void
copy_to(struct psc_mparticles *mprts, struct psc_mparticles *mprts_to,
	void (*put_particle)(struct vpic_mparticles_prt *prt, int n, void *ctx))
{
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts)->vmprts;
  
  int n_prts_by_patch[mprts->nr_patches];
  vpic_mparticles_get_size_all(vmprts, mprts->nr_patches, n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .mprts = mprts_to, .p = p };
    vpic_mparticles_get_grid_nx_dx(vmprts, ctx.im, ctx.dx);
    for (int d = 0; d < 3; d++) {
      ctx.im[d] += 2; // add ghost points
    }
    ctx.dVi = 1.f / (ctx.dx[0] * ctx.dx[1] * ctx.dx[2]);
    vpic_mparticles_get_particles(vmprts, n_prts, off, put_particle, &ctx);

    off += n_prts;
  }
}

// ======================================================================
// conversion to "single"

static void
put_particle_single(struct vpic_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = _ctx;
  particle_single_t *part = psc_mparticles_single_get_one(ctx->mprts, ctx->p, n);

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
  part->kind    = prt->kind;
  part->pxi     = prt->ux[0];
  part->pyi     = prt->ux[1];
  part->pzi     = prt->ux[2];
  part->qni_wni = ppsc->kinds[prt->kind].q * prt->w * ctx->dVi;
}

static void
psc_mparticles_vpic_copy_from_single(struct psc_mparticles *mprts,
				    struct psc_mparticles *mprts_single, unsigned int flags)
{
  assert(0);
}

static void
psc_mparticles_vpic_copy_to_single(struct psc_mparticles *mprts,
				  struct psc_mparticles *mprts_single, unsigned int flags)
{
  copy_to(mprts, mprts_single, put_particle_single);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_methods

static struct mrc_obj_method psc_mparticles_vpic_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_mparticles_vpic_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mparticles_vpic_copy_from_single),
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_vpic_setup

static void
psc_mparticles_vpic_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);

  sub->vmprts = vpic_mparticles_create();
  vpic_mparticles_ctor_from_simulation(sub->vmprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_get_size_all

static void
psc_mparticles_vpic_get_size_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts)->vmprts;

  vpic_mparticles_get_size_all(vmprts, mprts->nr_patches, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_reserve_all

static void
psc_mparticles_vpic_reserve_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  int cur_n_prts_by_patch[mprts->nr_patches];
  psc_mparticles_vpic_get_size_all(mprts, cur_n_prts_by_patch);
  for (int p = 0; p < mprts->nr_patches; p++) {
    if (n_prts_by_patch[p] != cur_n_prts_by_patch[p]) {
      mprintf("psc_mparticles_vpic_reserve_all: %d -> %d\n",
	      cur_n_prts_by_patch[p], cur_n_prts_by_patch[p]);
    }
    assert(n_prts_by_patch[p] == cur_n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_resize_all

static void
psc_mparticles_vpic_resize_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  int cur_n_prts_by_patch[mprts->nr_patches];
  psc_mparticles_vpic_get_size_all(mprts, cur_n_prts_by_patch);
  if (n_prts_by_patch[0] != cur_n_prts_by_patch[0]) {
    mprintf("psc_mparticles_vpic_resize_all: %d -> %d\n",
	    cur_n_prts_by_patch[0], n_prts_by_patch[0]);
  }
  assert(n_prts_by_patch[0] == cur_n_prts_by_patch[0]);
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_get_nr_particles

static unsigned int
psc_mparticles_vpic_get_nr_particles(struct psc_mparticles *mprts)
{
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts)->vmprts;
  
  return vpic_mparticles_get_nr_particles(vmprts);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "vpic"
  
struct psc_mparticles_ops psc_mparticles_vpic_ops = {
  .name                    = "vpic",
  .size                    = sizeof(struct psc_mparticles_vpic),
  .methods                 = psc_mparticles_vpic_methods,
  .setup                   = psc_mparticles_vpic_setup,
  .reserve_all             = psc_mparticles_vpic_reserve_all,
  .get_size_all            = psc_mparticles_vpic_get_size_all,
  .resize_all              = psc_mparticles_vpic_resize_all,
  .get_nr_particles        = psc_mparticles_vpic_get_nr_particles,
};

