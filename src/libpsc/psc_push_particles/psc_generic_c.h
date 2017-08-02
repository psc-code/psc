
#ifndef PSC_GENERIC_C_H
#define PSC_GENERIC_C_H

#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_c.h"

typedef particle_real_t creal;
#define creal_abs particle_real_abs
#define creal_sqrt particle_real_sqrt

void psc_push_particles_generic_c_push_a_y(struct psc_push_particles *push,
					   struct psc_particles *particles_base,
					   struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_z(struct psc_push_particles *push,
					   struct psc_particles *particles_base,
					   struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_xy(struct psc_push_particles *push,
					    struct psc_particles *particles_base,
					    struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_xz(struct psc_push_particles *push,
					    struct psc_particles *particles_base,
					    struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_yz(struct psc_push_particles *push,
					    struct psc_particles *particles_base,
					    struct psc_fields *flds_base);
void psc_push_particles_generic_c_push_a_xyz(struct psc_push_particles *push,
					     struct psc_particles *particles_base,
					     struct psc_fields *flds_base);

void psc_push_particles_generic_c_push_yz_a(struct psc_push_particles *push,
					    mparticles_base_t *particles_base,
					    mfields_base_t *flds_base);

#endif
