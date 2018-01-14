
#ifndef PSC_PARTICLES_AS_CUDA_H
#define PSC_PARTICLES_AS_CUDA_H

#include "psc_particles_cuda.h"

using mparticles_t = mparticles_cuda_t;
using particle_t = mparticles_t::particle_t;

#define PARTICLE_TYPE               "cuda"

#define PSC_PARTICLES_AS_CUDA 1

#include "particle_iter.h"

#endif

