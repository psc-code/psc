
#include "psc.h"
#include "psc_cuda.h"
#include "psc_particles_cuda.h"

#include <mrc_profile.h>

// FIXME -> header
void cuda_sort_patch(int p, particles_cuda_t *pp);

static void
particles_cuda_alloc(int p, particles_cuda_t *pp, int n_part,
		     bool need_block_offsets, bool need_cell_offsets)
{
  struct psc_patch *patch = &ppsc->patch[p];

  pp->n_part = n_part;
  int bs[3] = { BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z };
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      bs[d] = 1;
    }
    pp->b_mx[d] = (patch->ldims[d] + bs[d] - 1) / bs[d];
  }
  pp->nr_blocks = pp->b_mx[0] * pp->b_mx[1] * pp->b_mx[2];
  __particles_cuda_alloc(pp, need_block_offsets, need_cell_offsets);
}

static void
particles_cuda_free(particles_cuda_t *pp)
{
  __particles_cuda_free(pp);
}

#include "psc_particles_as_c.h"

static inline void
find_cell(real xi, real yi, real zi, int l[3])
{
  l[0] = cuda_nint(xi / ppsc->dx[0]);
  l[1] = cuda_nint(yi / ppsc->dx[1]);
  l[2] = cuda_nint(zi / ppsc->dx[2]);
  //  printf("l %d %d %d\n", l[0], l[1], l[2]);
}

// FIXME, should go away and always be done within cuda for consistency

static inline int
find_cellIdx(struct psc_patch *patch, struct cell_map *map,
	     particle_t *p)
{
  particle_real_t dxi = 1.f / ppsc->dx[0];
  particle_real_t dyi = 1.f / ppsc->dx[1];
  particle_real_t dzi = 1.f / ppsc->dx[2];
  particle_real_t xi[3] = {
    (p->xi - patch->xb[0]) * dxi,
    (p->yi - patch->xb[1]) * dyi,
    (p->zi - patch->xb[2]) * dzi };
  int pos[3];
  for (int d = 0; d < 3; d++) {
    pos[d] = cuda_fint(xi[d]);
  }
  
  return cell_map_3to1(map, pos);
}

static inline int
find_blockIdx(struct psc_patch *patch, struct cell_map *map,
	      particle_t *p)
{
  int cell_idx = find_cellIdx(patch, map, p);
  return cell_idx / (BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z);
}

static inline void
blockIdx_to_blockCrd(struct psc_patch *patch, struct cell_map *map,
		     int bidx, int bi[3])
{
  int cidx = bidx * (BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z);
  cell_map_1to3(map, cidx, bi);
  bi[0] /= BLOCKSIZE_X;
  bi[1] /= BLOCKSIZE_Y;
  bi[2] /= BLOCKSIZE_Z;
}

static void
psc_mparticles_copy_cf_to_cuda(mparticles_cuda_t *particles, mparticles_t *particles_cf,
			       bool need_block_offsets, bool need_cell_offsets)
{
  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = &ppsc->patch[p];
    particles_t *pp_cf = psc_mparticles_get_patch(particles_cf, p);
    particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
    pp->n_part = pp_cf->n_part;

    float4 *xi4  = calloc(pp_cf->n_part, sizeof(float4));
    float4 *pxi4 = calloc(pp_cf->n_part, sizeof(float4));

    for (int n = 0; n < pp_cf->n_part; n++) {
      particle_t *part_cf = particles_get_one(pp_cf, n);

      real qni = part_cf->qni;
      real wni = part_cf->wni;
      real qni_div_mni = qni / part_cf->mni;
      real qni_wni;
      if (qni != 0.) {
	qni_wni = qni * wni;
      } else {
	qni_wni = wni;
      }
      
      xi4[n].x  = part_cf->xi - patch->xb[0];
      xi4[n].y  = part_cf->yi - patch->xb[1];
      xi4[n].z  = part_cf->zi - patch->xb[2];
      xi4[n].w  = qni_div_mni;
      pxi4[n].x = part_cf->pxi;
      pxi4[n].y = part_cf->pyi;
      pxi4[n].z = part_cf->pzi;
      pxi4[n].w = qni_wni;
    }

    int bs[3] = { BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z };
    for (int d = 0; d < 3; d++) {
      if (ppsc->domain.gdims[d] == 1) {
	bs[d] = 1;
      }
      if (bs[d] != 1) {
	bs[d] *= 2; // sort not only within blocks, but also on lowest block
	// bit, so we can do the checkerboard passes
      }
    }
    struct cell_map map;
    cell_map_init(&map, patch->ldims, bs);

    int *offsets = NULL;
    if (0 && need_block_offsets) {
      // FIXME, should go away and can be taken over by c_offsets
      offsets = calloc(pp->nr_blocks + 1, sizeof(*offsets));
      int last_block = -1;
      for (int n = 0; n <= pp->n_part; n++) {
	int block;
	if (n < pp->n_part) {
	  particle_t *part_cf = particles_get_one(pp_cf, n);
	  block = find_blockIdx(patch, &map, part_cf);
	} else {
	  block = pp->nr_blocks;
	}
	assert(last_block <= block);
	while (last_block < block) {
	  offsets[last_block+1] = n;
	  last_block++;
	}
      }

#if 0
      for (int b = 0; b < pp->nr_blocks; b++) {
	int bi[3];
	blockIdx_to_blockCrd(patch, &map, b, bi);
	printf("block %d [%d,%d,%d]: %d:%d\n", b, bi[0], bi[1], bi[2],
	       offsets[b], offsets[b+1]);
      }
#endif
    }

    // FIXME, could be computed on the cuda side
    int *c_pos = calloc(map.N * 3, sizeof(*c_pos));
    for (int cidx = 0; cidx < map.N; cidx++) {
      int ci[3];
      cell_map_1to3(&map, cidx, ci);
      c_pos[3*cidx + 0] = ci[0];
      c_pos[3*cidx + 1] = ci[1];
      c_pos[3*cidx + 2] = ci[2];
    }

    int *c_offsets = NULL;
    if (need_cell_offsets) {
      const int cells_per_block = BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z;
      c_offsets = calloc(pp->nr_blocks * cells_per_block + 1,
			      sizeof(*c_offsets));
      int last_block = -1;
      for (int n = 0; n <= pp->n_part; n++) {
	int block;
	if (n < pp->n_part) {
	  particle_t *part_cf = particles_get_one(pp_cf, n);
	  block = find_cellIdx(patch, &map, part_cf);
	} else {
	  block = map.N;
	}
	assert(block <= pp->nr_blocks * cells_per_block);
	assert(last_block <= block);
	while (last_block < block) {
	  c_offsets[last_block+1] = n;
	  last_block++;
	}
      }
    }
    cell_map_free(&map);

    __particles_cuda_to_device(pp, xi4, pxi4, offsets, c_offsets, c_pos);

    if (need_block_offsets) {
      cuda_sort_patch(p, pp);
    }

    free(offsets);
    free(c_pos);
    free(c_offsets);
    free(xi4);
    free(pxi4);
  }
}

static void
psc_mparticles_copy_cf_from_cuda(mparticles_cuda_t *particles, mparticles_t *particles_cf)
{
  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = &ppsc->patch[p];
    particles_t *pp_cf = psc_mparticles_get_patch(particles_cf, p);
    particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
    pp_cf->n_part = pp->n_part;

    float4 *xi4  = calloc(pp->n_part, sizeof(float4));
    float4 *pxi4 = calloc(pp->n_part, sizeof(float4));

    __particles_cuda_from_device(pp, xi4, pxi4);

    for (int n = 0; n < pp_cf->n_part; n++) {
      particle_real_t qni_div_mni = xi4[n].w;
      particle_real_t qni_wni = pxi4[n].w;
      particle_real_t qni, mni, wni;
      if (qni_div_mni == 0.) {
	qni = 0.;
	wni = qni_wni;
	mni = -1.;
	// FIXME, irrelevant if no-copy assert(0); // can't recover the mass of a neutral particle
      } else {
	qni = qni_div_mni > 0 ? 1. : -1.;
	mni = qni / qni_div_mni;
	wni = qni_wni / qni;
      }

      particle_t *part_base = particles_get_one(pp_cf, n);
      part_base->xi  = xi4[n].x + patch->xb[0];
      part_base->yi  = xi4[n].y + patch->xb[1];
      part_base->zi  = xi4[n].z + patch->xb[2];
      part_base->pxi = pxi4[n].x;
      part_base->pyi = pxi4[n].y;
      part_base->pzi = pxi4[n].z;
      part_base->qni = qni;
      part_base->mni = mni;
      part_base->wni = wni;
    }

    free(xi4);
    free(pxi4);
  }
}

mparticles_cuda_t *
psc_mparticles_cuda_get_cuda(void *_particles_base, unsigned int flags)
{
  return _particles_base;
}

void
psc_mparticles_cuda_put_cuda(mparticles_cuda_t *particles, void *particles_base)
{
}

static bool __gotten;

mparticles_c_t *
psc_mparticles_cuda_get_c(void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_get_c", 1., 0, 0);
  }
  prof_start(pr);

  assert(!__gotten);
  __gotten = true;
    
  mparticles_base_t *particles_base = _particles_base;
  int *nr_particles_by_patch = malloc(particles_base->nr_patches * sizeof(int));
  for (int p = 0; p < particles_base->nr_patches; p++) {
    nr_particles_by_patch[p] = psc_mparticles_base_nr_particles_by_patch(particles_base, p);
  }
  struct mrc_domain *domain = particles_base->domain;
  mparticles_c_t *particles_c = psc_mparticles_c_create(mrc_domain_comm(domain));
  psc_mparticles_c_set_domain_nr_particles(particles_c, domain, nr_particles_by_patch);
  psc_mparticles_c_setup(particles_c);
  free(nr_particles_by_patch);

  psc_mparticles_copy_cf_from_cuda(particles_base, particles_c);

  prof_stop(pr);

  return particles_c;
}

void
psc_mparticles_cuda_put_c(mparticles_c_t *particles_c, void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_put_c", 1., 0, 0);
  }
  prof_start(pr);

  assert(__gotten);
  __gotten = false;

  mparticles_base_t *particles_base = _particles_base;
  psc_mparticles_copy_cf_to_cuda(particles_base, particles_c,
				 true, false); // FIXME, need to sort
  psc_mparticles_c_destroy(particles_c);

  prof_stop(pr);
}

mparticles_fortran_t *
psc_mparticles_cuda_get_fortran(void *_particles_base)
{
  assert(0);
}

void
psc_mparticles_cuda_put_fortran(mparticles_fortran_t *particles, void *_particles_base)
{
  assert(0);
}

// ======================================================================

mparticles_cuda_t *
psc_mparticles_c_get_cuda(void *_particles_base, unsigned int flags)
{
  assert(ppsc->nr_patches == 1); // many things would break...

  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_get_cuda", 1., 0, 0);
  }
  prof_start(pr);

  assert(!__gotten);
  __gotten = true;
  
  mparticles_c_t *particles_base = _particles_base;

  int *nr_particles_by_patch = malloc(particles_base->nr_patches * sizeof(int));
  for (int p = 0; p < particles_base->nr_patches; p++) {
    nr_particles_by_patch[p] = psc_mparticles_c_nr_particles_by_patch(particles_base, p);
  }
  struct mrc_domain *domain = particles_base->domain;
  mparticles_cuda_t *particles = psc_mparticles_cuda_create(mrc_domain_comm(domain));
  psc_mparticles_cuda_set_domain_nr_particles(particles, domain, nr_particles_by_patch);
  psc_mparticles_cuda_setup(particles);
  free(nr_particles_by_patch);

  psc_mparticles_copy_cf_to_cuda(particles, particles_base,
				 flags & MP_NEED_BLOCK_OFFSETS, flags & MP_NEED_CELL_OFFSETS);

  prof_stop(pr);
  return particles;
}

void
psc_mparticles_c_put_cuda(mparticles_cuda_t *particles, void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_cuda_put", 1., 0, 0);
  }
  prof_start(pr);

  assert(__gotten);
  __gotten = false;

  mparticles_t *particles_base = _particles_base;
  psc_mparticles_copy_cf_from_cuda(particles, particles_base);
  psc_mparticles_cuda_destroy(particles);

  prof_stop(pr);
}

// ======================================================================
// psc_mparticles_cuda

static void
_psc_mparticles_cuda_set_domain_nr_particles(mparticles_cuda_t *mparticles,
					     struct mrc_domain *domain,
					     int *nr_particles_by_patch)
{
  mparticles->domain = domain;
  mrc_domain_get_patches(domain, &mparticles->nr_patches);

  mparticles->data = calloc(mparticles->nr_patches, sizeof(*mparticles->data));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_cuda_alloc(p, psc_mparticles_get_patch_cuda(mparticles, p),
			 nr_particles_by_patch[p],
			 true, false); // FIXME, don't hardcode
  }
}

static int
_psc_mparticles_cuda_nr_particles_by_patch(mparticles_cuda_t *mparticles, int p)
{
  return psc_mparticles_get_patch_cuda(mparticles, p)->n_part;
}

static void
_psc_mparticles_cuda_destroy(mparticles_cuda_t *mparticles)
{
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_cuda_free(psc_mparticles_get_patch_cuda(mparticles, p));
  }
  free(mparticles->data);
}

// ======================================================================
// psc_mparticles: subclass "cuda"
  
struct psc_mparticles_cuda_ops psc_mparticles_cuda_ops = {
  .name                    = "cuda",
  .set_domain_nr_particles = _psc_mparticles_cuda_set_domain_nr_particles,
  .nr_particles_by_patch   = _psc_mparticles_cuda_nr_particles_by_patch,
};

static void
psc_mparticles_cuda_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mparticles_cuda, &psc_mparticles_cuda_ops);
}

struct mrc_class_psc_mparticles_cuda mrc_class_psc_mparticles_cuda = {
  .name             = "psc_mparticles_cuda",
  .size             = sizeof(struct psc_mparticles_cuda),
  .init             = psc_mparticles_cuda_init,
  .destroy          = _psc_mparticles_cuda_destroy,
};


