
#include "psc_cuda.h"

void
psc_bnd_cuda_exchange_particles(struct psc *psc, mparticles_cuda_t *particles)
{
  assert(psc->nr_patches == 1);
  int size;
  MPI_Comm_size(psc_comm(psc), &size);
  assert(size == 1);
  
  cuda_exchange_particles(0, &particles->p[0]);
}

