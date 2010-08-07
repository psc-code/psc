
#include "psc_cuda.h"

#include <stdlib.h>
#include <assert.h>

static inline void
find_cell(real xi, real yi, real zi, int l[3])
{
  l[0] = nint(xi / psc.dx[0]);
  l[1] = nint(yi / psc.dx[1]);
  l[2] = nint(zi / psc.dx[2]);
  //  printf("l %d %d %d\n", l[0], l[1], l[2]);
}

static inline int
find_blockIdx(struct psc_cuda *cuda, real xi, real yi, real zi)
{
  int bi[3];
  find_cell(xi, yi, zi, bi);
  bi[0] /= BLOCKSIZE_X;
  bi[1] /= BLOCKSIZE_Y;
  bi[2] /= BLOCKSIZE_Z;

  assert(bi[0] >= 0 && bi[0] < cuda->b_mx[0]);
  assert(bi[1] >= 0 && bi[1] < cuda->b_mx[1]);
  assert(bi[2] >= 0 && bi[2] < cuda->b_mx[2]);

  return (bi[2] * cuda->b_mx[1] + bi[1]) * cuda->b_mx[0] + bi[0];
}

static inline void
blockIdx_to_blockCrd(struct psc_cuda *cuda, int bidx, int bi[3])
{
  bi[2] = bidx / (cuda->b_mx[1] * cuda->b_mx[0]);
  bidx -= bi[2] * (cuda->b_mx[1] * cuda->b_mx[0]);
  bi[1] = bidx / cuda->b_mx[0];
  bidx -= bi[1] * cuda->b_mx[0];
  bi[0] = bidx;
}

// ======================================================================

static void
cuda_create()
{
  struct psc_cuda *cuda = malloc(sizeof(*cuda));
  psc.c_ctx = cuda;
}

static void
cuda_destroy()
{
  struct psc_cuda *cuda = psc.c_ctx;
  free(cuda);
}

static void
cuda_particles_from_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;

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

  cuda->xi4 = xi4;
  cuda->pxi4 = pxi4;

  cuda->b_mx[0] = (psc.ihi[0] - psc.ilo[0] + BLOCKSIZE_X - 1) / BLOCKSIZE_X;
  cuda->b_mx[1] = (psc.ihi[1] - psc.ilo[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y;
  cuda->b_mx[2] = (psc.ihi[2] - psc.ilo[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  cuda->nr_blocks = cuda->b_mx[0] * cuda->b_mx[1] * cuda->b_mx[2];
  cuda->offsets = calloc(cuda->nr_blocks + 1, sizeof(*cuda->offsets));
  int last_block = -1;
  for (int i = 0; i <= psc.pp.n_part; i++) {
    int block;
    if (i < psc.pp.n_part) {
      block = find_blockIdx(cuda, psc.pp.particles[i].xi, psc.pp.particles[i].yi,
			    psc.pp.particles[i].zi);
    } else {
      block = cuda->nr_blocks;
    }
    assert(last_block <= block);
    while (last_block < block) {
      cuda->offsets[last_block+1] = i;
      last_block++;
    }
  }

#if 0
  for (int c = 0; c < cuda->nr_blocks; c++) {
    int ci[3];
    blockIdx_to_blockCrd(cuda, c, ci);
    printf("cell %d [%d,%d,%d]: %d:%d\n", c, ci[0], ci[1], ci[2],
	   cuda->offsets[c], cuda->offsets[c+1]);
  }
#endif

  __cuda_particles_from_fortran(cuda);

  free(cuda->offsets);
}

static void
cuda_particles_to_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;
  float4 *xi4  = cuda->xi4;
  float4 *pxi4 = cuda->pxi4;

  __cuda_particles_to_fortran(cuda);

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

static void
cuda_fields_from_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;
  fields_cuda_t *pf = &cuda->f;

  pf->flds = calloc(NR_FIELDS * psc.fld_size, sizeof(*pf->flds));
  
  for (int m = EX; m <= HZ; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_CUDA(pf, m, jx,jy,jz) = F3_BASE(m, jx,jy,jz);
	}
      }
    }
  }

  __cuda_fields_from_fortran(pf);
}

static void
cuda_fields_to_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;
  fields_cuda_t *pf = &cuda->f;

  __cuda_fields_to_fortran(pf);

  free(pf->flds);
}

struct psc_ops psc_ops_cuda = {
  .name = "cuda",
  .create                 = cuda_create,
  .destroy                = cuda_destroy,
  .particles_from_fortran = cuda_particles_from_fortran,
  .particles_to_fortran   = cuda_particles_to_fortran,
  .fields_from_fortran    = cuda_fields_from_fortran,
  .fields_to_fortran      = cuda_fields_to_fortran,
  .push_part_yz_a         = cuda_push_part_yz_a,
  .push_part_yz_b         = cuda_push_part_yz_b2,
};
