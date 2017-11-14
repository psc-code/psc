
#include "vpic_mparticles.h"

#include <cassert>

extern vpic_simulation *simulation;

// ======================================================================
// vpic_mparticles

// ----------------------------------------------------------------------
// vpic_mparticles_create

struct vpic_mparticles *
vpic_mparticles_create()
{
  return new vpic_mparticles;
}

// ----------------------------------------------------------------------
// vpic_mparticles_ctor_from_simulation

void
vpic_mparticles_ctor_from_simulation(struct vpic_mparticles *vmprts)
{
  vmprts->species_list = simulation->species_list;
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_nr_particles

int vpic_mparticles_get_nr_particles(struct vpic_mparticles *vmprts)
{
  int n_prts = 0;

  species_t *sp;
  LIST_FOR_EACH(sp, vmprts->species_list) {
    n_prts += sp->np;
  }

  return n_prts;
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_size_all

void vpic_mparticles_get_size_all(struct vpic_mparticles *vmprts, int n_patches,
				  int *n_prts_by_patch)
{
  assert(n_patches == 1);
  n_prts_by_patch[0] = vpic_mparticles_get_nr_particles(vmprts);
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_grid_nx_dx

void vpic_mparticles_get_grid_nx_dx(struct vpic_mparticles *vmprts, int *nx, float *dx)
{
  nx[0] = vmprts->species_list->g->nx;
  nx[1] = vmprts->species_list->g->ny;
  nx[2] = vmprts->species_list->g->nz;
  dx[0] = vmprts->species_list->g->dx;
  dx[1] = vmprts->species_list->g->dy;
  dx[2] = vmprts->species_list->g->dz;
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_particles

void vpic_mparticles_get_particles(struct vpic_mparticles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*put_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx)
{
  species_t *sp;
  unsigned int v_off = 0;
  LIST_FOR_EACH(sp, vmprts->species_list) {
    unsigned int v_n_prts = sp->np;

    unsigned int nb = std::max(v_off, off), ne = std::min(v_off + v_n_prts, off + n_prts);
    for (int n = nb; n < ne; n++) {
      particle *p = &sp->p[n - v_off];
      struct vpic_mparticles_prt prt;
      prt.dx[0] = p->dx;
      prt.dx[1] = p->dy;
      prt.dx[2] = p->dz;
      prt.i     = p->i;
      prt.ux[0] = p->ux;
      prt.ux[1] = p->uy;
      prt.ux[2] = p->uz;
      prt.w     = p->w;
      prt.kind  = sp->id;
      put_particle(&prt, n - off, ctx);
    }

    v_off += v_n_prts;
  }
}

void vpic_mparticles_set_particles(struct vpic_mparticles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*get_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx)
{
  species_t *sp;
  unsigned int v_off = 0;
  LIST_FOR_EACH(sp, vmprts->species_list) {
    unsigned int v_n_prts = sp->np;

    unsigned int nb = std::max(v_off, off), ne = std::min(v_off + v_n_prts, off + n_prts);
    for (int n = nb; n < ne; n++) {
      struct vpic_mparticles_prt prt;
      get_particle(&prt, n - off, ctx);
      particle *p = &sp->p[n - v_off];
      p->dx = prt.dx[0];
      p->dy = prt.dx[1];
      p->dz = prt.dx[2];
      p->i  = prt.i;
      p->ux = prt.ux[0];
      p->uy = prt.ux[1];
      p->uz = prt.ux[2];
      p->w  = prt.w;
      assert(prt.kind == sp->id);
    }

    v_off += v_n_prts;
  }
}

// ----------------------------------------------------------------------
// conversions

struct copy_ctx {
  bk_mparticles *bkmprts;
  int p;
};

static void
copy_from(vpic_mparticles *vmprts, bk_mparticles *bkmprts,
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
copy_to(vpic_mparticles *vmprts, bk_mparticles *bkmprts,
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
get_particle_single_by_kind(struct vpic_mparticles_prt *prt, int n, void *_ctx)
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
put_particle_single_by_kind(struct vpic_mparticles_prt *prt, int n, void *_ctx)
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

void vpic_mparticles_copy_from_single_by_kind(vpic_mparticles *vmprts, bk_mparticles *bkmprts)
{
  copy_from(vmprts, bkmprts, get_particle_single_by_kind);
}

void vpic_mparticles_copy_to_single_by_kind(vpic_mparticles *vmprts, bk_mparticles *bkmprts)
{
  copy_to(vmprts, bkmprts, put_particle_single_by_kind);
}


