
#include "psc_cuda.h"

#include <stdlib.h>
#include <assert.h>

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

  float4 *xi4  = calloc(psc.n_part, sizeof(float4));
  float4 *pxi4 = calloc(psc.n_part, sizeof(float4));

  for (int i = 0; i < psc.n_part; i++) {
    real qni = psc.f_part[i].qni;
    real wni = psc.f_part[i].wni;
    real qni_div_mni = qni / psc.f_part[i].mni;
    real qni_wni;
    if (qni != 0.) {
      qni_wni = qni * wni;
    } else {
      qni_wni = wni;
    }

    xi4[i].x  = psc.f_part[i].xi;
    xi4[i].y  = psc.f_part[i].yi;
    xi4[i].z  = psc.f_part[i].zi;
    xi4[i].w  = qni_div_mni;
    pxi4[i].x = psc.f_part[i].pxi;
    pxi4[i].y = psc.f_part[i].pyi;
    pxi4[i].z = psc.f_part[i].pzi;
    pxi4[i].w = qni_wni;
  }

  cuda->xi4 = xi4;
  cuda->pxi4 = pxi4;

  __cuda_particles_from_fortran(cuda);
}

static void
cuda_particles_to_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;
  float4 *xi4  = cuda->xi4;
  float4 *pxi4 = cuda->pxi4;

  __cuda_particles_to_fortran(cuda);

  for (int i = 0; i < psc.n_part; i++) {
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

    psc.f_part[i].xi  = xi4[i].x;
    psc.f_part[i].yi  = xi4[i].y;
    psc.f_part[i].zi  = xi4[i].z;
    psc.f_part[i].pxi = pxi4[i].x;
    psc.f_part[i].pyi = pxi4[i].y;
    psc.f_part[i].pzi = pxi4[i].z;
    psc.f_part[i].qni = qni;
    psc.f_part[i].mni = mni;
    psc.f_part[i].wni = wni;
  }

  free(xi4);
  free(pxi4);
}

static void
cuda_fields_from_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;

  cuda->flds = calloc(NR_FIELDS * psc.fld_size, sizeof(float));
  
  for (int m = EX; m <= BZ; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  CF3(m, jx,jy,jz) = FF3(m, jx,jy,jz);
	}
      }
    }
  }

  __cuda_fields_from_fortran(cuda);
}

static void
cuda_fields_to_fortran()
{
  struct psc_cuda *cuda = psc.c_ctx;

  __cuda_fields_to_fortran(cuda);

  free(cuda->flds);
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
  .push_part_yz_b         = cuda_push_part_yz_b,
};
