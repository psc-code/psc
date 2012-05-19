
#include "psc.h"
#include "psc_cuda.h"
#include "psc_bnd_cuda.h"
#include "psc_particles_cuda.h"
#include "psc_push_particles.h"

EXTERN_C void cuda_init(int rank);

static void *
_psc_mparticles_cuda_alloc_patch(int p, int n_part, unsigned int flags)
{
  if (p == 0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cuda_init(rank);
  }

  struct psc_patch *patch = &ppsc->patch[p];
  particles_cuda_t *pp = calloc(1, sizeof(*pp));

  pp->n_part = n_part;
  int bs[3];
  for (int d = 0; d < 3; d++) {
    switch (flags & MP_BLOCKSIZE_MASK) {
    case MP_BLOCKSIZE_1X1X1: bs[d] = 1; break;
    case MP_BLOCKSIZE_2X2X2: bs[d] = 2; break;
    case MP_BLOCKSIZE_4X4X4: bs[d] = 4; break;
    case MP_BLOCKSIZE_8X8X8: bs[d] = 8; break;
    default: assert(0);
    }
    if (ppsc->domain.gdims[d] == 1) {
      bs[d] = 1;
    }
    pp->blocksize[d] = bs[d];
    assert(patch->ldims[d] % bs[d] == 0); // not sure what breaks if not
    pp->b_mx[d] = (patch->ldims[d] + bs[d] - 1) / bs[d];
    pp->b_dxi[d] = 1.f / (pp->blocksize[d] * ppsc->dx[d]);
  }
  pp->nr_blocks = pp->b_mx[0] * pp->b_mx[1] * pp->b_mx[2];

  for (int d = 0; d < 3; d++) {
    if (flags & MP_NO_CHECKERBOARD) {
      bs[d] = 1;
    } else {
      bs[d] = (patch->ldims[d] == 1) ? 1 : 2;
    }
  }
  cell_map_init(&pp->map, pp->b_mx, bs);

  __particles_cuda_alloc(pp, true, true); // FIXME, need separate flags

  cuda_alloc_block_indices(pp, &pp->d_part.bidx); // FIXME, merge into ^^^
  cuda_alloc_block_indices(pp, &pp->d_part.ids);
  cuda_alloc_block_indices(pp, &pp->d_part.alt_bidx);
  cuda_alloc_block_indices(pp, &pp->d_part.alt_ids);
  cuda_alloc_block_indices(pp, &pp->d_part.sums);

  pp->d_part.sort_ctx = sort_pairs_create(pp->b_mx);

  return pp;
}

static void
_psc_mparticles_cuda_free_patch(mparticles_cuda_t *mparticles, int p)
{
  particles_cuda_t *pp = psc_mparticles_get_patch_cuda(mparticles, p);

  cuda_free_block_indices(pp->d_part.bidx);
  cuda_free_block_indices(pp->d_part.ids);
  cuda_free_block_indices(pp->d_part.alt_bidx);
  cuda_free_block_indices(pp->d_part.alt_ids);
  cuda_free_block_indices(pp->d_part.sums);
  sort_pairs_destroy(pp->d_part.sort_ctx);
  __particles_cuda_free(pp);
  cell_map_free(&pp->map);
}

#include "psc_particles_as_c.h"

// FIXME, should go away and always be done within cuda for consistency

static inline int
find_cellIdx(struct psc_patch *patch, struct cell_map *map,
	     particles_t *pp, int n)
{
  particle_t *p = particles_get_one(pp, n);
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
	      particles_t *pp, int n, int blocksize[3])
{
  int cell_idx = find_cellIdx(patch, map, pp, n);
  return cell_idx / (blocksize[0] * blocksize[1] * blocksize[2]);
}

static inline void
blockIdx_to_blockCrd(struct psc_patch *patch, struct cell_map *map,
		     int bidx, int bi[3], int blocksize[3])
{
  int cidx = bidx * (blocksize[0] * blocksize[1] * blocksize[2]);
  cell_map_1to3(map, cidx, bi);
  for (int d = 0; d < 3; d++) {
    bi[d] /= blocksize[d];
  }
}

static void
_psc_mparticles_cuda_copy_from_c(int p, mparticles_cuda_t *particles,
				 mparticles_t *particles_cf, unsigned int flags)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particles_t *pp_cf = psc_mparticles_get_patch(particles_cf, p);
  particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
  pp->n_part = pp_cf->n_part;
  assert(pp->n_part <= pp->n_alloced);
  
  float4 *xi4  = calloc(pp_cf->n_part, sizeof(float4));
  float4 *pxi4 = calloc(pp_cf->n_part, sizeof(float4));
  
  for (int n = 0; n < pp_cf->n_part; n++) {
    particle_t *part_cf = particles_get_one(pp_cf, n);
    
    real qni = part_cf->qni;
    real wni = part_cf->wni;
    // FIXME!!! KH hack
    if (part_cf->kind == 1 || part_cf->kind == 3) {
      wni *= (1 + 1e-6);
    }
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
    
    float xi[3] = { xi4[n].x, xi4[n].y, xi4[n].z };
    for (int d = 0; d < 3; d++) {
      int bi = cuda_fint(xi[d] * pp->b_dxi[d]);
      if (bi < 0 || bi >= pp->b_mx[d]) {
	printf("XXX p %d xi %g %g %g\n", p, xi[0], xi[1], xi[2]);
	printf("XXX p %d n %d d %d xi4[n] %g biy %d // %d\n",
	       p, n, d, xi[d], bi, pp->b_mx[d]);
	if (bi < 0) {
	  xi[d] = 0.f;
	} else {
	  xi[d] *= (1. - 1e-6);
	}
      }
      bi = cuda_fint(xi[d] * pp->b_dxi[d]);
      assert(bi >= 0 && bi < pp->b_mx[d]);
    }
    xi4[n].x = xi[0];
    xi4[n].y = xi[1];
    xi4[n].z = xi[2];
  }
  
  int bs[3];
  for (int d = 0; d < 3; d++) {
    bs[d] = pp->blocksize[d];
    if (bs[d] != 1) {
      bs[d] *= 2; // sort not only within blocks, but also on lowest block
      // bit, so we can do the checkerboard passes
    }
  }
  struct cell_map map;
  cell_map_init(&map, patch->ldims, bs); // FIXME, already have it elsewhere
  
  int *offsets = NULL;
  if (0 && (flags & MP_NEED_BLOCK_OFFSETS)) {
    // FIXME, should go away and can be taken over by c_offsets
    offsets = calloc(pp->nr_blocks + 1, sizeof(*offsets));
    int last_block = -1;
    for (int n = 0; n <= pp->n_part; n++) {
      int block;
      if (n < pp->n_part) {
	block = find_blockIdx(patch, &map, pp_cf, n, pp->blocksize);
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
      blockIdx_to_blockCrd(patch, &map, b, bi, pp->blocksize);
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
  if (0 && (flags & MP_NEED_CELL_OFFSETS)) {
    const int cells_per_block = pp->blocksize[0] * pp->blocksize[1] * pp->blocksize[2];
    c_offsets = calloc(pp->nr_blocks * cells_per_block + 1,
		       sizeof(*c_offsets));
    int last_block = -1;
    for (int n = 0; n <= pp->n_part; n++) {
      int block;
      if (n < pp->n_part) {
	block = find_cellIdx(patch, &map, pp_cf, n);
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
  
  if (particles->flags & MP_NEED_BLOCK_OFFSETS) {
    cuda_sort_patch(p, pp);
  }
  if (particles->flags & MP_NEED_CELL_OFFSETS) {
    cuda_sort_patch_by_cell(p, pp);
  }
  // FIXME, sorting twice because we need both would be suboptimal
  if ((particles->flags & MP_NEED_CELL_OFFSETS) && 
      (particles->flags & MP_NEED_BLOCK_OFFSETS)) {
    MHERE;
  }
  
  free(offsets);
  free(c_pos);
  free(c_offsets);
  free(xi4);
  free(pxi4);
}

static void
_psc_mparticles_cuda_copy_to_c(int p, mparticles_cuda_t *particles,
			       mparticles_c_t *particles_cf, unsigned int flags)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particles_t *pp_cf = psc_mparticles_get_patch(particles_cf, p);
  particles_cuda_t *pp = psc_mparticles_get_patch_cuda(particles, p);
  pp_cf->n_part = pp->n_part;
  assert(pp_cf->n_part <= pp_cf->n_alloced);
  
  float4 *xi4  = calloc(pp->n_part, sizeof(float4));
  float4 *pxi4 = calloc(pp->n_part, sizeof(float4));
  
  __particles_cuda_from_device(pp, xi4, pxi4);
  
  for (int n = 0; n < pp_cf->n_part; n++) {
    particle_real_t qni_div_mni = xi4[n].w;
    particle_real_t qni_wni = pxi4[n].w;
    particle_real_t qni, mni, wni;
    unsigned int kind;
    if (qni_div_mni == 0.) {
      qni = 0.;
      wni = qni_wni;
      mni = -1.;
      assert(0); // can't recover the mass of a neutral particle
    } else {
      qni = qni_div_mni > 0 ? 1. : -1.;
      mni = qni / qni_div_mni;
      wni = qni_wni / qni;
      kind = (qni > 0.) ? 2 : 0;
      if (wni > 1. + .5e-7) {
	kind++;
      }
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
    part_base->kind = kind;
  }

  free(xi4);
  free(pxi4);
}

// ======================================================================
// psc_mparticles_cuda

static int
_psc_mparticles_cuda_nr_particles_by_patch(mparticles_cuda_t *mparticles, int p)
{
  return psc_mparticles_get_patch_cuda(mparticles, p)->n_part;
}

// ======================================================================
// psc_mparticles: subclass "cuda"
  
static struct mrc_obj_method _psc_mparticles_cuda_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",         _psc_mparticles_cuda_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c",       _psc_mparticles_cuda_copy_from_c),
  {}
};

struct psc_mparticles_ops psc_mparticles_cuda_ops = {
  .name                    = "cuda",
  .methods                 = _psc_mparticles_cuda_methods,
  .nr_particles_by_patch   = _psc_mparticles_cuda_nr_particles_by_patch,
  .alloc_patch             = _psc_mparticles_cuda_alloc_patch,
  .free_patch              = _psc_mparticles_cuda_free_patch,
  .size_of_particles_t     = sizeof(particles_cuda_t),
};

