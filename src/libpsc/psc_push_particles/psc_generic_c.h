
#ifndef PSC_GENERIC_C_H
#define PSC_GENERIC_C_H

#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_c.h"

typedef particle_real_t creal;
#define creal_abs particle_real_abs
#define creal_sqrt particle_real_sqrt

void psc_push_particles_generic_c_push_mprts_y(struct psc_push_particles *push,
					       struct psc_mparticles *mprts,
					       struct psc_mfields *mflds);
void psc_push_particles_generic_c_push_mprts_z(struct psc_push_particles *push,
					       struct psc_mparticles *mprts,
					       struct psc_mfields *mflds);

void psc_push_particles_generic_c_push_mprts_xy(struct psc_push_particles *push,
						struct psc_mparticles *mprts,
						struct psc_mfields *mflds);
void psc_push_particles_generic_c_push_mprts_xz(struct psc_push_particles *push,
						struct psc_mparticles *mprts,
						struct psc_mfields *mflds);
void psc_push_particles_generic_c_push_mprts_yz(struct psc_push_particles *push,
						struct psc_mparticles *mprts,
						struct psc_mfields *mflds);

void psc_push_particles_generic_c_push_mprts_xyz(struct psc_push_particles *push,
						 struct psc_mparticles *mprts,
						 struct psc_mfields *mflds);

#endif
