
#ifndef PSC_PARTICLES_AS_FORTRAN_H
#define PSC_PARTICLES_AS_FORTRAN_H

#include "psc_particles_fortran.h"

typedef particle_fortran_real_t particle_real_t;
typedef particle_fortran_t particle_t;

#define mparticles_get_one            psc_mparticles_fortran_get_one
#define mparticles_get_n_prts         psc_mparticles_fortran_get_n_prts
#define mparticles_patch_reserve      psc_mparticles_fortran_patch_reserve
#define mparticles_patch_push_back    psc_mparticles_fortran_patch_push_back
#define mparticles_patch_resize       psc_mparticles_fortran_patch_resize
#define mparticles_patch_capacity     psc_mparticles_fortran_patch_capacity
#define mparticles_patch_get_b_mx     psc_mparticles_fortran_patch_get_b_mx
#define mparticles_patch_get_b_dxi    psc_mparticles_fortran_patch_get_b_dxi

#define particle_buf_t              psc_particle_fortran_buf_t
#define particle_buf_t              psc_particle_fortran_buf_t
#define particle_buf_ctor           psc_particle_fortran_buf_ctor
#define particle_buf_dtor           psc_particle_fortran_buf_dtor
#define particle_buf_size           psc_particle_fortran_buf_size
#define particle_buf_resize         psc_particle_fortran_buf_resize
#define particle_buf_reserve        psc_particle_fortran_buf_reserve
#define particle_buf_push_back      psc_particle_fortran_buf_push_back

#define particle_real_fint            particle_fortran_real_fint
#define particle_iter_t               psc_particle_fortran_iter_t
#define particle_iter_equal           psc_particle_fortran_iter_equal
#define particle_iter_next            psc_particle_fortran_iter_next 
#define particle_iter_deref           psc_particle_fortran_iter_deref 
#define particle_iter_at              psc_particle_fortran_iter_at
#define particle_range_t              psc_particle_fortran_range_t
#define particle_range_mprts          psc_particle_fortran_range_mprts
#define particle_range_size           psc_particle_fortran_range_size

#define MPI_PARTICLES_REAL            MPI_PARTICLES_FORTRAN_REAL
#define PARTICLE_TYPE                 "fortran"

#define PSC_PARTICLES_AS_FORTRAN 1

#include "particle_iter.h"

#endif

