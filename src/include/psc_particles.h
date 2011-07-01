
#ifndef PSC_PARTICLES_H
#define PSC_PARTICLES_H

// ----------------------------------------------------------------------
// base particles type

#if PARTICLES_BASE == PARTICLES_FORTRAN

typedef particles_fortran_t particles_base_t;
typedef mparticles_fortran_t mparticles_base_t;
typedef particle_fortran_t particle_base_t;
typedef particle_fortran_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL MPI_PARTICLES_FORTRAN_REAL

#define particles_base_alloc   particles_fortran_alloc
#define particles_base_realloc particles_fortran_realloc
#define particles_base_free    particles_fortran_free
#define particles_base_get_one particles_fortran_get_one
#define mparticles_base_create  mparticles_fortran_create
#define mparticles_base_set_domain_nr_particles mparticles_fortran_set_domain_nr_particles
#define mparticles_base_setup   mparticles_fortran_setup
#define mparticles_base_destroy mparticles_fortran_destroy

#elif PARTICLES_BASE == PARTICLES_C

#include "psc_particles_c.h"

typedef particles_c_t particles_base_t;
typedef mparticles_c_t mparticles_base_t;
typedef particle_c_t particle_base_t;
typedef particle_c_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL    MPI_PARTICLES_C_REAL

#define particles_base_alloc   particles_c_alloc
#define particles_base_realloc particles_c_realloc
#define particles_base_free    particles_c_free
#define particles_base_get_one particles_c_get_one
#define mparticles_base_create  mparticles_c_create
#define mparticles_base_set_domain_nr_particles mparticles_c_set_domain_nr_particles
#define mparticles_base_setup   mparticles_c_setup
#define mparticles_base_destroy mparticles_c_destroy

#elif PARTICLES_BASE == PARTICLES_SSE2

#include "psc_particles_sse2.h"

typedef particles_sse2_t particles_base_t;
typedef mparticles_sse2_t mparticles_base_t;
typedef particle_sse2_t particle_base_t;
typedef particle_sse2_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL    MPI_PARTICLES_SSE2_REAL

#define particles_base_alloc   particles_sse2_alloc
#define particles_base_realloc particles_sse2_realloc
#define particles_base_free    particles_sse2_free
#define particles_base_get_one particles_sse2_get_one
#define mparticles_base_create  mparticles_sse2_create
#define mparticles_base_set_domain_nr_particles mparticles_sse2_set_domain_nr_particles
#define mparticles_base_setup   mparticles_sse2_setup
#define mparticles_base_destroy mparticles_sse2_destroy

#elif PARTICLES_BASE == PARTICLES_CBE

#include "psc_particles_cbe.h"

typedef particles_cbe_t particles_base_t;
typedef mparticles_cbe_t mparticles_base_t;
typedef particle_cbe_t particle_base_t;
typedef particle_cbe_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL    MPI_PARTICLES_CBE_REAL

#define particles_base_alloc   particles_cbe_alloc
#define particles_base_realloc particles_cbe_realloc
#define particles_base_free    particles_cbe_free
#define particles_base_get_one particles_cbe_get_one
#define mparticles_base_create  mparticles_cbe_create
#define mparticles_base_set_domain_nr_particles mparticles_cbe_set_domain_nr_particles
#define mparticles_base_setup   mparticles_cbe_setup
#define mparticles_base_destroy mparticles_cbe_destroy

#else
#error unknown PARTICLES_BASE
#endif

mparticles_base_t *mparticles_base_create(MPI_Comm comm);
void mparticles_base_set_domain_nr_particles(mparticles_base_t *mparticles, 
					     struct mrc_domain *domain,
					     int *nr_particles_by_patch);
void mparticles_base_setup(mparticles_base_t *mparticles);
void mparticles_base_destroy(mparticles_base_t *mparticles);


#endif
