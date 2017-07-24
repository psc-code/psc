
#include "cuda_iface.h"
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <cstdio>

#define sqr(x) ((x)*(x))

static cuda_heating_foil foil;

// ----------------------------------------------------------------------
// cuda_heating_setup_foil

void
cuda_heating_setup_foil(struct cuda_heating_foil *_foil)
{
  foil = *_foil;

  float width = foil.zh - foil.zl;
  foil.fac = (8.f * pow(foil.T, 1.5)) / (sqrt(foil.Mi) * width);
}

// ----------------------------------------------------------------------
// foil_get_H

static float
foil_get_H(float *xx)
{
  if (xx[2] <= foil.zl || xx[2] >= foil.zh) {
    return 0;
  }

  return foil.fac * exp(-(sqr(xx[0] - foil.xc) +
			  sqr(xx[1] - foil.yc)) / sqr(foil.rH));
}

// ----------------------------------------------------------------------
// particle_kick

static float4
particle_kick(float4 pxi4, float H)
{
  float ran1, ran2, ran3, ran4, ran5, ran6;
  do {
    ran1 = random() / ((float) RAND_MAX + 1);
    ran2 = random() / ((float) RAND_MAX + 1);
    ran3 = random() / ((float) RAND_MAX + 1);
    ran4 = random() / ((float) RAND_MAX + 1);
    ran5 = random() / ((float) RAND_MAX + 1);
    ran6 = random() / ((float) RAND_MAX + 1);
  } while (ran1 >= 1.f || ran2 >= 1.f || ran3 >= 1.f ||
	   ran4 >= 1.f || ran5 >= 1.f || ran6 >= 1.f);

  float ranx = sqrtf(-2.f*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
  float rany = sqrtf(-2.f*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4);
  float ranz = sqrtf(-2.f*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6);

  float Dp = sqrtf(H * foil.heating_dt);

  pxi4.x += Dp * ranx;
  pxi4.y += Dp * rany;
  pxi4.z += Dp * ranz;

  return pxi4;
}

// ----------------------------------------------------------------------
// cuda_heating_run_foil_gold

void
cuda_heating_run_foil_gold(struct cuda_mparticles *cmprts)
{
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cmprts->d_pxi4);
  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);
  thrust::device_ptr<unsigned int> d_off(cmprts->d_off);

  for (int b = 0; b < cmprts->n_blocks; b++) {
    int p = b / cmprts->n_blocks_per_patch;
    for (int n = d_off[b]; n < d_off[b+1]; n++) {
      float4 xi4 = d_xi4[n];

      int prt_kind = cuda_float_as_int(xi4.w);
      if (prt_kind != foil.kind) {
	continue;
      }

      float *xb = &cmprts->xb_by_patch[p][0];
      float xx[3] = {
	xi4.x + xb[0],
	xi4.y + xb[1],
	xi4.z + xb[2],
      };

      float H = foil_get_H(xx);
      float4 pxi4 = d_pxi4[n];
      printf("%s xx = %g %g %g H = %g px = %g %g %g\n", (H > 0) ? "H" : " ",
	     xx[0], xx[1], xx[2], H,
	     pxi4.x, pxi4.y, pxi4.z);
      if (H > 0) {
	d_pxi4[n] = particle_kick(d_pxi4[n], H);
	// pxi4 = d_pxi4[n];
	// printf("H xx = %g %g %g H = %g px = %g %g %g\n", xx[0], xx[1], xx[2], H,
	//        pxi4.x, pxi4.y, pxi4.z);
      }
    }
  }
}

// ----------------------------------------------------------------------
// cuda_heating_run_foil

void
cuda_heating_run_foil(struct cuda_mparticles *cmprts)
{
  printf("cuda_heating_run_foil\n");
  cuda_heating_run_foil_gold(cmprts);
}

