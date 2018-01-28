
#include "cuda_particles_bnd_iface.h"
#include "cuda_mparticles.h"

void cuda_particles_bnd::prep(cuda_mparticles* cmprts)
{
  cmprts->bnd_prep();
}


