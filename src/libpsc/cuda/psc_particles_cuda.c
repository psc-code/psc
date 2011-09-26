
#include "psc.h"
#include "psc_particles_cuda.h"

#include <mrc_profile.h>

static inline void
find_cell(real xi, real yi, real zi, int l[3])
{
  l[0] = cuda_nint(xi / ppsc->dx[0]);
  l[1] = cuda_nint(yi / ppsc->dx[1]);
  l[2] = cuda_nint(zi / ppsc->dx[2]);
  //  printf("l %d %d %d\n", l[0], l[1], l[2]);
}

static inline int
find_cellIdx(struct psc_patch *patch, struct cell_map *map,
	     particle_base_t *p)
{
  particle_base_real_t dxi = 1.f / ppsc->dx[0];
  particle_base_real_t dyi = 1.f / ppsc->dx[1];
  particle_base_real_t dzi = 1.f / ppsc->dx[2];
  particle_base_real_t xi[3] = {
    (p->xi - patch->xb[0]) * dxi,
    (p->yi - patch->xb[1]) * dyi,
    (p->zi - patch->xb[2]) * dzi };
  int pos[3];
  for (int d = 0; d < 3; d++) {
    pos[d] = particle_base_real_nint(xi[d]);
  }
  
  return cell_map_3to1(map, pos);
}

static inline int
find_blockIdx(struct psc_patch *patch, struct cell_map *map,
	      particle_base_t *p)
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

// ======================================================================

static bool __gotten;

void
psc_mparticles_cuda_get_from(mparticles_cuda_t *particles, void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_cuda_get", 1., 0, 0);
  }
  prof_start(pr);

  assert(!__gotten);
  __gotten = true;
    
  mparticles_base_t *particles_base = _particles_base;

  particles->p = calloc(ppsc->nr_patches, sizeof(*particles->p));
  assert(ppsc->nr_patches == 1); // many things would break...
  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = &ppsc->patch[p];
    particles_base_t *pp_base = &particles_base->p[p];
    particles_cuda_t *pp = &particles->p[p];
    particles_cuda_dev_t *h_part = &pp->h_part;

    pp->n_part = pp_base->n_part;
    float4 *xi4  = calloc(pp->n_part, sizeof(float4));
    float4 *pxi4 = calloc(pp->n_part, sizeof(float4));

    for (int n = 0; n < pp_base->n_part; n++) {
      particle_base_t *part_base = particles_base_get_one(pp_base, n);

      real qni = part_base->qni;
      real wni = part_base->wni;
      real qni_div_mni = qni / part_base->mni;
      real qni_wni;
      if (qni != 0.) {
	qni_wni = qni * wni;
      } else {
	qni_wni = wni;
      }
      
      xi4[n].x  = part_base->xi - patch->xb[0];
      xi4[n].y  = part_base->yi - patch->xb[1];
      xi4[n].z  = part_base->zi - patch->xb[2];
      xi4[n].w  = qni_div_mni;
      pxi4[n].x = part_base->pxi;
      pxi4[n].y = part_base->pyi;
      pxi4[n].z = part_base->pzi;
      pxi4[n].w = qni_wni;
    }

    h_part->xi4 = xi4;
    h_part->pxi4 = pxi4;

    int bs[3] = { BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z };
    for (int d = 0; d < 3; d++) {
      if (ppsc->domain.gdims[d] == 1) {
	bs[d] = 1;
      }
      pp->b_mx[d] = (patch->ldims[d] + bs[d] - 1) / bs[d];
    }
    pp->nr_blocks = pp->b_mx[0] * pp->b_mx[1] * pp->b_mx[2];
    // FIXME, should go away and can be taken over by c_offsets
    h_part->offsets = calloc(pp->nr_blocks + 1, sizeof(*h_part->offsets));

    for (int d = 0; d < 3; d++) {
      if (bs[d] != 1) {
	bs[d] *= 2; // sort not only within blocks, but also on lowest block
	// bit, so we can do the checkerboard passes
      }
    }
    int last_block = -1;
    struct cell_map map;
    cell_map_init(&map, patch->ldims, bs);
    for (int n = 0; n <= pp->n_part; n++) {
      int block;
      if (n < pp->n_part) {
	particle_base_t *part_base = particles_base_get_one(pp_base, n);
	block = find_blockIdx(patch, &map, part_base);
      } else {
	block = pp->nr_blocks;
      }
      if (last_block > block) {
	MHERE;
	break;
      }
      assert(last_block <= block);
      while (last_block < block) {
	h_part->offsets[last_block+1] = n;
	last_block++;
      }
    }

#if 0
    MHERE;
    for (int b = 0; b < pp->nr_blocks; b++) {
      int bi[3];
      blockIdx_to_blockCrd(patch, &map, b, bi);
      printf("block %d [%d,%d,%d]: %d:%d\n", b, bi[0], bi[1], bi[2],
	     h_part->offsets[b], h_part->offsets[b+1]);
    }
#endif
    h_part->c_pos = calloc(map.N * 3, sizeof(int));
    for (int cidx = 0; cidx < map.N; cidx++) {
      int ci[3];
      cell_map_1to3(&map, cidx, ci);
      h_part->c_pos[3*cidx + 0] = ci[0];
      h_part->c_pos[3*cidx + 1] = ci[1];
      h_part->c_pos[3*cidx + 2] = ci[2];
    }

    const int cells_per_block = BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z;
    h_part->c_offsets = calloc(pp->nr_blocks * cells_per_block + 1,
			       sizeof(*h_part->c_offsets));
    last_block = -1;
    for (int n = 0; n <= pp->n_part; n++) {
      int block;
      if (n < pp->n_part) {
	particle_base_t *part_base = particles_base_get_one(pp_base, n);
	block = find_cellIdx(patch, &map, part_base);
      } else {
	block = map.N;
      }
      assert(block <= pp->nr_blocks * cells_per_block);
      if (last_block > block) {
	MHERE;
	break;
      }
      assert(last_block <= block);
      while (last_block < block) {
	h_part->c_offsets[last_block+1] = n;
	last_block++;
      }
    }
    cell_map_free(&map);
    
    __particles_cuda_get(pp);
    
    free(h_part->offsets); // FIXME?!!!
    free(h_part->c_pos);
    free(h_part->c_offsets);
  }

  prof_stop(pr);
}

void
psc_mparticles_cuda_put_to(mparticles_cuda_t *particles, void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_cuda_put", 1., 0, 0);
  }
  prof_start(pr);

  assert(__gotten);
  __gotten = false;

  mparticles_base_t *particles_base = _particles_base;
  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = &ppsc->patch[p];
    particles_base_t *pp_base = &particles_base->p[p];
    particles_cuda_t *pp = &particles->p[p];
    assert(pp->n_part == pp_base->n_part);

    particles_cuda_dev_t *h_part = &pp->h_part;
    float4 *xi4  = h_part->xi4;
    float4 *pxi4 = h_part->pxi4;

    __particles_cuda_put(pp);

    for (int n = 0; n < pp_base->n_part; n++) {
      particle_base_real_t qni_div_mni = xi4[n].w;
      particle_base_real_t qni_wni = pxi4[n].w;
      particle_base_real_t qni, mni, wni;
      if (qni_div_mni == 0.) {
	qni = 0.;
	wni = qni_wni;
	mni = -1.;
	assert(0); // can't recover the mass of a neutral particle
      } else {
	qni = qni_div_mni > 0 ? 1. : -1.;
	mni = qni / qni_div_mni;
	wni = qni_wni / qni;
      }

      particle_base_t *part_base = particles_base_get_one(pp_base, n);
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
  free(particles->p);
  particles->p = NULL;

  prof_stop(pr);
}

