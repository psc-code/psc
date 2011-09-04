
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

static unsigned int ldims_bits[3];

static unsigned int
find_cell_indices_4x4x4(int p)
{
  struct psc_patch *patch = &ppsc->patch[p];

  unsigned int N = 1;
  for (int d = 0; d < 3; d++) {
    ldims_bits[d] = 0; // FIXME, need to handle invariant dirs generically
    while (patch->ldims[d] > (1 << ldims_bits[d])) {
      ldims_bits[d]++;
    }
    // need to be powers of 2 currently, could (should) be relaxed
    // to just require divisible by 4
    //    printf("ld %d %d\n", patch->ldims[d], ldims_bits[d]);
    assert(patch->ldims[d] == (1 << ldims_bits[d]));
    N *= (1 << ldims_bits[d]);
  }
  //  printf("N %d\n", N);
  return N;
}

static inline int
get_cell_index_4x4x4(int p, real xi, real yi, real zi)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_base_real_t dxi = 1.f / ppsc->dx[0];
  particle_base_real_t dyi = 1.f / ppsc->dx[1];
  particle_base_real_t dzi = 1.f / ppsc->dx[2];
  int *ldims = patch->ldims;
  
  particle_base_real_t u = (xi - patch->xb[0]) * dxi;
  particle_base_real_t v = (yi - patch->xb[1]) * dyi;
  particle_base_real_t w = (zi - patch->xb[2]) * dzi;
  unsigned int j0 = particle_base_real_nint(u);
  unsigned int j1 = particle_base_real_nint(v);
  unsigned int j2 = particle_base_real_nint(w);
  assert(j0 < ldims[0]);
  assert(j1 < ldims[1]);
  assert(j2 < ldims[2]);
  assert(ldims[0] == 1);

  return (((j2 >> 2) << (ldims_bits[1] - 2 + 4)) |
	  ((j1 >> 2) << (4)) |
	  (((j2 >> 1) & 1) << 3) |
	  (((j1 >> 1) & 1) << 2) |
	  ((j2 & 1) << 1) |
	  ((j1 & 1) << 0));
}

static inline int
find_cellIdx(particles_cuda_t *pp, real xi, real yi, real zi)
{
  int p = 0; // FIXME
  static bool first_time = true;
  if (first_time) {
    first_time = false;
    int N = find_cell_indices_4x4x4(p); // FIXME, once per patch
    assert(N == BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z * 
	   pp->b_mx[0] * pp->b_mx[1] * pp->b_mx[2]);
  }
  return get_cell_index_4x4x4(p, xi, yi, zi);
}

static inline int
find_blockIdx(particles_cuda_t *pp, real xi, real yi, real zi)
{
#if 0

  int cell_idx = find_cellIdx(pp, xi, yi, zi);
#if BLOCKSIZE_X == 1 && BLOCKSIZE_Y == 4 && BLOCKSIZE_Z == 4
  return cell_idx / 16;
#else
#error TBD
#endif

#else
  int bi[3];
  find_cell(xi, yi, zi, bi);
  bi[0] /= BLOCKSIZE_X;
  bi[1] /= BLOCKSIZE_Y;
  bi[2] /= BLOCKSIZE_Z;

  assert(bi[0] >= 0 && bi[0] < pp->b_mx[0]);
  assert(bi[1] >= 0 && bi[1] < pp->b_mx[1]);
  assert(bi[2] >= 0 && bi[2] < pp->b_mx[2]);

  return (bi[2] * pp->b_mx[1] + bi[1]) * pp->b_mx[0] + bi[0];
#endif
}

static inline void
blockIdx_to_blockCrd(particles_cuda_t *pp, int bidx, int bi[3])
{
  bi[2] = bidx / (pp->b_mx[1] * pp->b_mx[0]);
  bidx -= bi[2] * (pp->b_mx[1] * pp->b_mx[0]);
  bi[1] = bidx / pp->b_mx[0];
  bidx -= bi[1] * pp->b_mx[0];
  bi[0] = bidx;
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
      
      xi4[n].x  = part_base->xi;
      xi4[n].y  = part_base->yi;
      xi4[n].z  = part_base->zi;
      xi4[n].w  = qni_div_mni;
      pxi4[n].x = part_base->pxi;
      pxi4[n].y = part_base->pyi;
      pxi4[n].z = part_base->pzi;
      pxi4[n].w = qni_wni;
    }

    h_part->xi4 = xi4;
    h_part->pxi4 = pxi4;

    pp->b_mx[0] = (patch->ldims[0] + BLOCKSIZE_X - 1) / BLOCKSIZE_X;
    pp->b_mx[1] = (patch->ldims[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y;
    pp->b_mx[2] = (patch->ldims[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
    pp->nr_blocks = pp->b_mx[0] * pp->b_mx[1] * pp->b_mx[2];
    // FIXME, should go away and can be taken over by c_offsets
    h_part->offsets = calloc(pp->nr_blocks + 1, sizeof(*h_part->offsets));
    int last_block = -1;
    for (int n = 0; n <= pp->n_part; n++) {
      int block;
      if (n < pp->n_part) {
	particle_base_t *part_base = particles_base_get_one(pp_base, n);
	block = find_blockIdx(pp, part_base->xi, part_base->yi, part_base->zi);
      } else {
	block = pp->nr_blocks;
      }
      assert(last_block <= block);
      while (last_block < block) {
	h_part->offsets[last_block+1] = n;
	last_block++;
      }
    }

#if 0
    for (int c = 0; c < pp->nr_blocks; c++) {
      int ci[3];
      blockIdx_to_blockCrd(pp, c, ci);
      printf("cell %d [%d,%d,%d]: %d:%d\n", c, ci[0], ci[1], ci[2],
	     h_part->offsets[c], h_part->offsets[c+1]);
    }
#endif

    const int cells_per_block = BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z;
    h_part->c_offsets = calloc(pp->nr_blocks * cells_per_block + 1,
			       sizeof(*h_part->c_offsets));
    last_block = -1;
    for (int n = 0; n <= pp->n_part; n++) {
      int block;
      if (n < pp->n_part) {
	particle_base_t *part_base = particles_base_get_one(pp_base, n);
	block = find_cellIdx(pp, part_base->xi, part_base->yi, part_base->zi);
      } else {
	block = pp->nr_blocks * cells_per_block;
      }
      assert(block <= pp->nr_blocks * cells_per_block);
      assert(last_block <= block);
      while (last_block < block) {
	h_part->c_offsets[last_block+1] = n;
	last_block++;
      }
    }
    
    __particles_cuda_get(pp);
    
    free(h_part->offsets); // FIXME?!!!
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
      part_base->xi  = xi4[n].x;
      part_base->yi  = xi4[n].y;
      part_base->zi  = xi4[n].z;
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

