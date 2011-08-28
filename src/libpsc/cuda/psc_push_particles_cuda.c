
#include "psc_push_particles_private.h"
#include "psc_cuda.h"

// ======================================================================
// psc_push_particles: subclass "cuda"

struct psc_push_particles_ops psc_push_particles_cuda_ops = {
  .name                  = "cuda",
  .push_yz_a             = psc_push_particles_cuda_push_yz_a,
  .push_yz_b             = psc_push_particles_cuda_push_yz_b,
  .push_yz               = psc_push_particles_cuda_push_yz,
};
