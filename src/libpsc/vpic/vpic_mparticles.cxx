
#include "psc_particles.h"
#include <cassert>

#include "vpic_config.h"
#include "vpic_iface.h"

// ======================================================================
// vpic_mparticles

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

// ----------------------------------------------------------------------
// conversions

struct copy_ctx {
  bk_mparticles *bkmprts;
  int p;
};

static void
copy_from(Particles *vmprts, bk_mparticles *bkmprts,
	  void (*get_particle)(struct vpic_mparticles_prt *prt, int n, void *ctx))
{
  int n_patches = 1; // FIXME
  int n_prts_by_patch[n_patches];
  vpic_mparticles_get_size_all(vmprts, n_patches, n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .bkmprts = bkmprts, .p = p };
    vpic_mparticles_set_particles(vmprts, n_prts, off, get_particle, &ctx);

    off += n_prts;
  }
}
 
static void
copy_to(Particles *vmprts, bk_mparticles *bkmprts,
	void (*put_particle)(struct vpic_mparticles_prt *prt, int n, void *ctx))
{
  int n_patches = 1; // FIXME
  int n_prts_by_patch[n_patches];
  vpic_mparticles_get_size_all(vmprts, n_patches, n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .bkmprts = bkmprts, .p = p };
    vpic_mparticles_get_particles(vmprts, n_prts, off, put_particle, &ctx);

    off += n_prts;
  }
}

static void
convert_from_single_by_kind(struct vpic_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = (struct copy_ctx *) _ctx;
  particle_single_by_kind *part = &ctx->bkmprts->at(ctx->p, n);

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
  struct copy_ctx *ctx = (struct copy_ctx *) _ctx;
  particle_single_by_kind *part = &ctx->bkmprts->at(ctx->p, n);

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

void vpic_mparticles_copy_from_single_by_kind(Particles *vmprts, bk_mparticles *bkmprts)
{
  copy_from(vmprts, bkmprts, convert_from_single_by_kind);
}

void vpic_mparticles_copy_to_single_by_kind(Particles *vmprts, bk_mparticles *bkmprts)
{
  copy_to(vmprts, bkmprts, convert_to_single_by_kind);
}


