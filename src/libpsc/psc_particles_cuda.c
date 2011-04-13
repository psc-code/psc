
#include "psc.h"
#include "psc_particles_cuda.h"

static inline void
find_cell(real xi, real yi, real zi, int l[3])
{
  l[0] = nint(xi / psc.dx[0]);
  l[1] = nint(yi / psc.dx[1]);
  l[2] = nint(zi / psc.dx[2]);
  //  printf("l %d %d %d\n", l[0], l[1], l[2]);
}

static inline int
find_blockIdx(particles_cuda_t *pp, real xi, real yi, real zi)
{
  int bi[3];
  find_cell(xi, yi, zi, bi);
  bi[0] /= BLOCKSIZE_X;
  bi[1] /= BLOCKSIZE_Y;
  bi[2] /= BLOCKSIZE_Z;

  assert(bi[0] >= 0 && bi[0] < pp->b_mx[0]);
  assert(bi[1] >= 0 && bi[1] < pp->b_mx[1]);
  assert(bi[2] >= 0 && bi[2] < pp->b_mx[2]);

  return (bi[2] * pp->b_mx[1] + bi[1]) * pp->b_mx[0] + bi[0];
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

void
particles_cuda_get(particles_cuda_t *pp)
{
  particles_cuda_dev_t *h_part = &pp->h_part;

  float4 *xi4  = calloc(psc.pp.n_part, sizeof(float4));
  float4 *pxi4 = calloc(psc.pp.n_part, sizeof(float4));

  for (int i = 0; i < psc.pp.n_part; i++) {
    //    particle_base_t *p = particles_base_get_one(&psc.pp, i);
    real qni = psc.pp.particles[i].qni;
    real wni = psc.pp.particles[i].wni;
    real qni_div_mni = qni / psc.pp.particles[i].mni;
    real qni_wni;
    if (qni != 0.) {
      qni_wni = qni * wni;
    } else {
      qni_wni = wni;
    }

    xi4[i].x  = psc.pp.particles[i].xi;
    xi4[i].y  = psc.pp.particles[i].yi;
    xi4[i].z  = psc.pp.particles[i].zi;
    xi4[i].w  = qni_div_mni;
    pxi4[i].x = psc.pp.particles[i].pxi;
    pxi4[i].y = psc.pp.particles[i].pyi;
    pxi4[i].z = psc.pp.particles[i].pzi;
    pxi4[i].w = qni_wni;
  }

  h_part->xi4 = xi4;
  h_part->pxi4 = pxi4;

  pp->b_mx[0] = (psc.ihi[0] - psc.ilo[0] + BLOCKSIZE_X - 1) / BLOCKSIZE_X;
  pp->b_mx[1] = (psc.ihi[1] - psc.ilo[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y;
  pp->b_mx[2] = (psc.ihi[2] - psc.ilo[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  pp->nr_blocks = pp->b_mx[0] * pp->b_mx[1] * pp->b_mx[2];
  h_part->offsets = calloc(pp->nr_blocks + 1, sizeof(*h_part->offsets));
  int last_block = -1;
  for (int i = 0; i <= psc.pp.n_part; i++) {
    int block;
    if (i < psc.pp.n_part) {
      block = find_blockIdx(pp, psc.pp.particles[i].xi, psc.pp.particles[i].yi,
			    psc.pp.particles[i].zi);
    } else {
      block = pp->nr_blocks;
    }
    assert(last_block <= block);
    while (last_block < block) {
      h_part->offsets[last_block+1] = i;
      last_block++;
    }
  }

#if 0
  for (int c = 0; c < cuda->nr_blocks; c++) {
    int ci[3];
    blockIdx_to_blockCrd(cuda, c, ci);
    printf("cell %d [%d,%d,%d]: %d:%d\n", c, ci[0], ci[1], ci[2],
	   h_part->offsets[c], h_part->offsets[c+1]);
  }
#endif

  __particles_cuda_get(pp);

  free(h_part->offsets);
}

void
particles_cuda_put(particles_cuda_t *pp)
{
  particles_cuda_dev_t *h_part = &pp->h_part;
  float4 *xi4  = h_part->xi4;
  float4 *pxi4 = h_part->pxi4;

  __particles_cuda_put(pp);

  for (int i = 0; i < psc.pp.n_part; i++) {
    f_real qni_div_mni = xi4[i].w;
    f_real qni_wni = pxi4[i].w;
    f_real qni, mni, wni;
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

    psc.pp.particles[i].xi  = xi4[i].x;
    psc.pp.particles[i].yi  = xi4[i].y;
    psc.pp.particles[i].zi  = xi4[i].z;
    psc.pp.particles[i].pxi = pxi4[i].x;
    psc.pp.particles[i].pyi = pxi4[i].y;
    psc.pp.particles[i].pzi = pxi4[i].z;
    psc.pp.particles[i].qni = qni;
    psc.pp.particles[i].mni = mni;
    psc.pp.particles[i].wni = wni;
  }

  free(xi4);
  free(pxi4);
}

