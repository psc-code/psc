
#include "psc_particles.h"
#include <cassert>

#include "vpic_iface.h"

// ======================================================================
// vpic_mparticles

// ----------------------------------------------------------------------
// vpic_mparticles_new_from_simulation

Particles *
vpic_mparticles_new_from_simulation(Simulation *sim)
{
  return &sim->particles_;
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_nr_particles

int vpic_mparticles_get_nr_particles(Particles *vmprts)
{
  int n_prts = 0;

  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    n_prts += sp->np;
  }
  // species_t *sp;
  // LIST_FOR_EACH(sp, vmprts->sl_) {
  //   n_prts += sp->np;
  // }

  return n_prts;
}

// ----------------------------------------------------------------------
// vpic_mparticles_reserve_all
//
// This is a bit iffy, since we don't really want to reallocate stuff here,
// at least for now, and we wouldn't be able to know how to split this into
// the different species, anyway.

void vpic_mparticles_reserve_all(Particles *vmprts, int n_patches,
				 int *n_prts_by_patch)
{
  assert(n_patches == 1);

  for (int p = 0; p < n_patches; p++) {
    int n_prts = 0, n_prts_alloced = 0;
    for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
      n_prts += sp->np;
      n_prts_alloced += sp->max_np;
    }
#if 0
    if (n_prts_by_patch[p] != n_prts) {
      mprintf("vpic_mparticles_reserve_all: %d (currently %d max %d)\n",
	      n_prts_by_patch[p], n_prts, n_prts_alloced);
    }
#endif
    assert(n_prts_by_patch[p] <= n_prts_alloced);
  }
}

// ----------------------------------------------------------------------
// vpic_mparticles_resize_all
//
// Even more iffy, since can't really resize the per-species arrays, since we don't
// know how the total # of particles we're given should be divided up

void vpic_mparticles_resize_all(Particles *vmprts, int n_patches,
				int *n_prts_by_patch)
{
  assert(n_patches == 1);
  
  // we can't resize to the numbers given, unless it's "resize to 0", we'll just do nothing
  // The mparticles conversion function should call resize_all() itself first, resizing to
  // 0, and then using push_back, which will increase the count back to the right value

  if (n_prts_by_patch[0] == 0) {
    for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
      sp->np = 0;
    }
  } else {
#if 0
    int cur_n_prts_by_patch[n_patches];
    vpic_mparticles_get_size_all(vmprts, n_patches, cur_n_prts_by_patch);

    mprintf("vpic_mparticles_resize_all: ignoring %d -> %d\n",
	    cur_n_prts_by_patch[0], n_prts_by_patch[0]);
#endif
  }
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_size_all

void vpic_mparticles_get_size_all(Particles *vmprts, int n_patches,
				  int *n_prts_by_patch)
{
  assert(n_patches == 1);
  n_prts_by_patch[0] = vpic_mparticles_get_nr_particles(vmprts);
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_grid_nx_dx

void vpic_mparticles_get_grid_nx_dx(Particles *vmprts, int *nx, float *dx)
{
  grid_t *g = vmprts->getGrid_t();
  nx[0] = g->nx;
  nx[1] = g->ny;
  nx[2] = g->nz;
  dx[0] = g->dx;
  dx[1] = g->dy;
  dx[2] = g->dz;
}

// ----------------------------------------------------------------------
// vpic_mparticles_get_particles

void vpic_mparticles_get_particles(Particles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*put_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx)
{
  unsigned int v_off = 0;
  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    unsigned int v_n_prts = sp->np;

    unsigned int nb = std::max(v_off, off), ne = std::min(v_off + v_n_prts, off + n_prts);
    for (int n = nb; n < ne; n++) {
      particle *p = &sp->p[n - v_off];
#if 0
      int i = p->i;
      int im[3] = { sp->g->nx + 2, sp->g->ny + 2, sp->g->nz + 2 };
      int i3[3];
      i3[2] = i / (im[0] * im[1]); i -= i3[2] * (im[0] * im[1]);
      i3[1] = i / im[0]; i-= i3[1] * im[0];
      i3[0] = i;
      if (!(i3[2] >= 1 && i3[2] <= sp->g->nz)) {
	mprintf("i3 %d %d %d\n", i3[0], i3[1], i3[2]);
	assert(0);
      }
#endif
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

void vpic_mparticles_set_particles(Particles *vmprts, unsigned int n_prts, unsigned int off,
				   void (*get_particle)(struct vpic_mparticles_prt *, int, void *),
				   void *ctx)
{
  unsigned int v_off = 0;
  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
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

void vpic_mparticles_push_back(Particles *vmprts, const struct vpic_mparticles_prt *prt)
{
  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    if (sp->id == prt->kind) {
      assert(sp->np < sp->max_np);
      // the below is inject_particle_raw()
      particle_t * RESTRICT p = sp->p + (sp->np++);
      p->dx = prt->dx[0]; p->dy = prt->dx[1]; p->dz = prt->dx[2]; p->i = prt->i;
      p->ux = prt->ux[0]; p->uy = prt->ux[1]; p->uz = prt->ux[2]; p->w = prt->w;
      return;
    }
  }
  mprintf("prt->kind %d not found in species list!\n", prt->kind);
  assert(0);
}

// ----------------------------------------------------------------------
// vpic_mparticles_sort

void vpic_mparticles_sort(Particles *vmprts, int step)
{
  // Sort the particles for performance if desired.
  
  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    if (sp->sort_interval > 0 && (step % sp->sort_interval) == 0) {
      mpi_printf(MPI_COMM_WORLD, "Performance sorting \"%s\"\n", sp->name);
      TIC vmprts->sort_p(&*sp); TOC(sort_p, 1);
    }
  }
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

void vpic_mparticles_copy_from_single_by_kind(Particles *vmprts, bk_mparticles *bkmprts)
{
  copy_from(vmprts, bkmprts, get_particle_single_by_kind);
}

void vpic_mparticles_copy_to_single_by_kind(Particles *vmprts, bk_mparticles *bkmprts)
{
  copy_to(vmprts, bkmprts, put_particle_single_by_kind);
}


