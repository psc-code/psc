
#include "psc.h"
#include "psc_cuda.h"
#include "psc_bnd_cuda.h"
#include "psc_particles_cuda.h"
#include "psc_particles_single.h"
#include "psc_push_particles.h"

EXTERN_C void cuda_init(int rank);

// ======================================================================
// psc_particles "cuda"

static void
psc_particles_cuda_setup(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

  if (prts->p == 0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cuda_init(rank);
  }

  struct psc_patch *patch = &ppsc->patch[prts->p];

  int bs[3];
  for (int d = 0; d < 3; d++) {
    switch (prts->flags & MP_BLOCKSIZE_MASK) {
    case MP_BLOCKSIZE_1X1X1: bs[d] = 1; break;
    case MP_BLOCKSIZE_2X2X2: bs[d] = 2; break;
    case MP_BLOCKSIZE_4X4X4: bs[d] = 4; break;
    case MP_BLOCKSIZE_8X8X8: bs[d] = 8; break;
    default: assert(0);
    }
    if (ppsc->domain.gdims[d] == 1) {
      bs[d] = 1;
    }
    cuda->blocksize[d] = bs[d];
    assert(patch->ldims[d] % bs[d] == 0); // not sure what breaks if not
    cuda->b_mx[d] = (patch->ldims[d] + bs[d] - 1) / bs[d];
    cuda->b_dxi[d] = 1.f / (cuda->blocksize[d] * ppsc->dx[d]);
  }
  cuda->nr_blocks = cuda->b_mx[0] * cuda->b_mx[1] * cuda->b_mx[2];

  for (int d = 0; d < 3; d++) {
    if (prts->flags & MP_NO_CHECKERBOARD) {
      bs[d] = 1;
    } else {
      bs[d] = (patch->ldims[d] == 1) ? 1 : 2;
    }
  }
  cell_map_init(&cuda->map, cuda->b_mx, bs);

  __particles_cuda_alloc(prts, true, true); // FIXME, need separate flags

  cuda_alloc_block_indices(prts, &cuda->d_part.bidx); // FIXME, merge into ^^^
  cuda_alloc_block_indices(prts, &cuda->d_part.ids);
  cuda_alloc_block_indices(prts, &cuda->d_part.alt_bidx);
  cuda_alloc_block_indices(prts, &cuda->d_part.alt_ids);
  cuda_alloc_block_indices(prts, &cuda->d_part.sums);

  cuda->d_part.sort_ctx = sort_pairs_create(cuda->b_mx);
}

static void
psc_particles_cuda_destroy(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

  cuda_free_block_indices(cuda->d_part.bidx);
  cuda_free_block_indices(cuda->d_part.ids);
  cuda_free_block_indices(cuda->d_part.alt_bidx);
  cuda_free_block_indices(cuda->d_part.alt_ids);
  cuda_free_block_indices(cuda->d_part.sums);
  sort_pairs_destroy(cuda->d_part.sort_ctx);
  __particles_cuda_free(prts);
  cell_map_free(&cuda->map);
}

// FIXME, should go away and always be done within cuda for consistency

static inline int
find_cellIdx(struct psc_patch *patch, struct cell_map *map,
	     struct psc_particles *pp, int n)
{
  particle_c_t *p = particles_c_get_one(pp, n);
  particle_c_real_t dxi = 1.f / ppsc->dx[0];
  particle_c_real_t dyi = 1.f / ppsc->dx[1];
  particle_c_real_t dzi = 1.f / ppsc->dx[2];
  particle_c_real_t xi[3] = { p->xi * dxi, p->yi * dyi, p->zi * dzi };
  int pos[3];
  for (int d = 0; d < 3; d++) {
    pos[d] = cuda_fint(xi[d]);
  }
  
  return cell_map_3to1(map, pos);
}

static inline int
find_blockIdx(struct psc_patch *patch, struct cell_map *map,
	      struct psc_particles *pp, int n, int blocksize[3])
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

// ======================================================================
// conversion to "c"

static void
psc_particles_cuda_copy_from_c(struct psc_particles *prts_cuda,
			       struct psc_particles *prts_c, unsigned int flags)
{
  int p = prts_cuda->p;
  struct psc_patch *patch = &ppsc->patch[p];
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts_cuda);
  prts_cuda->n_part = prts_c->n_part;
  assert(prts_cuda->n_part <= cuda->n_alloced);
  
  float4 *xi4  = calloc(prts_c->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts_c->n_part, sizeof(float4));
  
  for (int n = 0; n < prts_c->n_part; n++) {
    particle_c_t *part_c = particles_c_get_one(prts_c, n);
    
    real qni = part_c->qni;
    real wni = part_c->wni;
    // FIXME!!! KH hack
    if (part_c->kind == 1 || part_c->kind == 3) {
      wni *= (1 + 1e-6);
    }
    real qni_div_mni = qni / part_c->mni;
    real qni_wni;
    if (qni != 0.) {
      qni_wni = qni * wni;
    } else {
      qni_wni = wni;
    }
    
    xi4[n].x  = part_c->xi;
    xi4[n].y  = part_c->yi;
    xi4[n].z  = part_c->zi;
    xi4[n].w  = qni_div_mni;
    pxi4[n].x = part_c->pxi;
    pxi4[n].y = part_c->pyi;
    pxi4[n].z = part_c->pzi;
    pxi4[n].w = qni_wni;
    
    float xi[3] = { xi4[n].x, xi4[n].y, xi4[n].z };
    for (int d = 0; d < 3; d++) {
      int bi = cuda_fint(xi[d] * cuda->b_dxi[d]);
      if (bi < 0 || bi >= cuda->b_mx[d]) {
	printf("XXX p %d xi %g %g %g\n", p, xi[0], xi[1], xi[2]);
	printf("XXX p %d n %d d %d xi4[n] %g biy %d // %d\n",
	       p, n, d, xi[d], bi, cuda->b_mx[d]);
	if (bi < 0) {
	  xi[d] = 0.f;
	} else {
	  xi[d] *= (1. - 1e-6);
	}
      }
      bi = cuda_fint(xi[d] * cuda->b_dxi[d]);
      assert(bi >= 0 && bi < cuda->b_mx[d]);
    }
    xi4[n].x = xi[0];
    xi4[n].y = xi[1];
    xi4[n].z = xi[2];
  }
  
  int bs[3];
  for (int d = 0; d < 3; d++) {
    bs[d] = cuda->blocksize[d];
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
    offsets = calloc(cuda->nr_blocks + 1, sizeof(*offsets));
    int last_block = -1;
    for (int n = 0; n <= prts_cuda->n_part; n++) {
      int block;
      if (n < prts_cuda->n_part) {
	block = find_blockIdx(patch, &map, prts_c, n, cuda->blocksize);
      } else {
	block = cuda->nr_blocks;
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
    const int cells_per_block = cuda->blocksize[0] * cuda->blocksize[1] * cuda->blocksize[2];
    c_offsets = calloc(cuda->nr_blocks * cells_per_block + 1, sizeof(*c_offsets));
    int last_block = -1;
    for (int n = 0; n <= prts_cuda->n_part; n++) {
      int block;
      if (n < prts_cuda->n_part) {
	block = find_cellIdx(patch, &map, prts_c, n);
      } else {
	block = map.N;
      }
      assert(block <= cuda->nr_blocks * cells_per_block);
      assert(last_block <= block);
      while (last_block < block) {
	c_offsets[last_block+1] = n;
	last_block++;
      }
    }
  }
  cell_map_free(&map);
  
  __particles_cuda_to_device(prts_cuda, xi4, pxi4, offsets, c_offsets, c_pos);
  
  if (prts_cuda->flags & MP_NEED_BLOCK_OFFSETS) {
    cuda_sort_patch(p, prts_cuda);
  }
  if (prts_cuda->flags & MP_NEED_CELL_OFFSETS) {
    cuda_sort_patch_by_cell(p, prts_cuda);
  }
  // FIXME, sorting twice because we need both would be suboptimal
  if ((prts_cuda->flags & MP_NEED_CELL_OFFSETS) && 
      (prts_cuda->flags & MP_NEED_BLOCK_OFFSETS)) {
    MHERE;
  }
  
  free(offsets);
  free(c_pos);
  free(c_offsets);
  free(xi4);
  free(pxi4);
}

static void
psc_particles_cuda_copy_to_c(struct psc_particles *prts_cuda,
			     struct psc_particles *prts_c, unsigned int flags)
{
  struct psc_particles_c *c = psc_particles_c(prts_c);
  prts_c->n_part = prts_cuda->n_part;
  assert(prts_c->n_part <= c->n_alloced);
  
  float4 *xi4  = calloc(prts_cuda->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts_cuda->n_part, sizeof(float4));
  
  __particles_cuda_from_device(prts_cuda, xi4, pxi4);
  
  for (int n = 0; n < prts_c->n_part; n++) {
    particle_c_real_t qni_div_mni = xi4[n].w;
    particle_c_real_t qni_wni = pxi4[n].w;
    particle_c_real_t qni, mni, wni;
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
    
    particle_c_t *part_base = particles_c_get_one(prts_c, n);
    part_base->xi  = xi4[n].x;
    part_base->yi  = xi4[n].y;
    part_base->zi  = xi4[n].z;
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
// conversion to "single"

static void
psc_particles_cuda_copy_from_single(struct psc_particles *prts_cuda,
				    struct psc_particles *prts, unsigned int flags)
{
  int p = prts_cuda->p;
  struct psc_patch *patch = &ppsc->patch[p];
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts_cuda);
  prts_cuda->n_part = prts->n_part;
  assert(prts_cuda->n_part <= cuda->n_alloced);
  
  float4 *xi4  = calloc(prts->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts->n_part, sizeof(float4));
  
  for (int n = 0; n < prts->n_part; n++) {
    particle_single_t *part = particles_single_get_one(prts, n);
    
    real qni_wni = part->qni_wni;
    // FIXME!!! KH hack
    if (part->kind == 1 || part->kind == 3) {
      qni_wni *= (1 + 1e-6);
    }
    
    xi4[n].x  = part->xi;
    xi4[n].y  = part->yi;
    xi4[n].z  = part->zi;
    xi4[n].w  = ppsc->kinds[part->kind].q / ppsc->kinds[part->kind].m;
    pxi4[n].x = part->pxi;
    pxi4[n].y = part->pyi;
    pxi4[n].z = part->pzi;
    pxi4[n].w = qni_wni;
    
    float xi[3] = { xi4[n].x, xi4[n].y, xi4[n].z };
    for (int d = 0; d < 3; d++) {
      int bi = cuda_fint(xi[d] * cuda->b_dxi[d]);
      if (bi < 0 || bi >= cuda->b_mx[d]) {
	printf("XXX p %d xi %g %g %g\n", p, xi[0], xi[1], xi[2]);
	printf("XXX p %d n %d d %d xi4[n] %g biy %d // %d\n",
	       p, n, d, xi[d], bi, cuda->b_mx[d]);
	if (bi < 0) {
	  xi[d] = 0.f;
	} else {
	  xi[d] *= (1. - 1e-6);
	}
      }
      bi = cuda_fint(xi[d] * cuda->b_dxi[d]);
      assert(bi >= 0 && bi < cuda->b_mx[d]);
    }
    xi4[n].x = xi[0];
    xi4[n].y = xi[1];
    xi4[n].z = xi[2];
  }
  
  int bs[3];
  for (int d = 0; d < 3; d++) {
    bs[d] = cuda->blocksize[d];
    if (bs[d] != 1) {
      bs[d] *= 2; // sort not only within blocks, but also on lowest block
      // bit, so we can do the checkerboard passes
    }
  }
  struct cell_map map;
  cell_map_init(&map, patch->ldims, bs); // FIXME, already have it elsewhere
  
  // FIXME, could be computed on the cuda side
  int *c_pos = calloc(map.N * 3, sizeof(*c_pos));
  for (int cidx = 0; cidx < map.N; cidx++) {
    int ci[3];
    cell_map_1to3(&map, cidx, ci);
    c_pos[3*cidx + 0] = ci[0];
    c_pos[3*cidx + 1] = ci[1];
    c_pos[3*cidx + 2] = ci[2];
  }
  
  cell_map_free(&map);
  
  __particles_cuda_to_device(prts_cuda, xi4, pxi4, NULL, NULL, c_pos);
  
  if (prts_cuda->flags & MP_NEED_BLOCK_OFFSETS) {
    cuda_sort_patch(p, prts_cuda);
  }
  if (prts_cuda->flags & MP_NEED_CELL_OFFSETS) {
    cuda_sort_patch_by_cell(p, prts_cuda);
  }
  // FIXME, sorting twice because we need both would be suboptimal
  if ((prts_cuda->flags & MP_NEED_CELL_OFFSETS) && 
      (prts_cuda->flags & MP_NEED_BLOCK_OFFSETS)) {
    MHERE;
  }
  
  free(c_pos);
  free(xi4);
  free(pxi4);
}

static void
psc_particles_cuda_copy_to_single(struct psc_particles *prts_cuda,
				  struct psc_particles *prts, unsigned int flags)
{
  struct psc_particles_single *sngl = psc_particles_single(prts);
  prts->n_part = prts_cuda->n_part;
  assert(prts->n_part <= sngl->n_alloced);
  
  float4 *xi4  = calloc(prts_cuda->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts_cuda->n_part, sizeof(float4));
  
  __particles_cuda_from_device(prts_cuda, xi4, pxi4);
  
  for (int n = 0; n < prts->n_part; n++) {
    particle_c_real_t qni_div_mni = xi4[n].w;
    particle_c_real_t qni_wni = pxi4[n].w;
    particle_c_real_t qni, mni, wni;
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
      wni = 1.f;
    }
    
    particle_single_t *part_base = particles_single_get_one(prts, n);
    part_base->xi  = xi4[n].x;
    part_base->yi  = xi4[n].y;
    part_base->zi  = xi4[n].z;
    part_base->pxi = pxi4[n].x;
    part_base->pyi = pxi4[n].y;
    part_base->pzi = pxi4[n].z;
    part_base->qni_wni = qni * wni;
    part_base->kind = kind;
  }

  free(xi4);
  free(pxi4);
}

// ======================================================================
// psc_mparticles: subclass "cuda"
  
struct psc_mparticles_ops psc_mparticles_cuda_ops = {
  .name                    = "cuda",
};

// ======================================================================
// psc_particles: subclass "cuda"

static struct mrc_obj_method psc_particles_cuda_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_particles_cuda_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_particles_cuda_copy_from_c),
  MRC_OBJ_METHOD("copy_to_single"  , psc_particles_cuda_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_particles_cuda_copy_from_single),
  {}
};

struct psc_particles_ops psc_particles_cuda_ops = {
  .name                    = "cuda",
  .size                    = sizeof(struct psc_particles_cuda),
  .methods                 = psc_particles_cuda_methods,
  .setup                   = psc_particles_cuda_setup,
  .destroy                 = psc_particles_cuda_destroy,
};
