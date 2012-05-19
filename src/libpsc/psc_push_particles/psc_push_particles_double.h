
#ifndef PSC_PUSH_PARTICLES_DOUBLE_H
#define PSC_PUSH_PARTICLES_DOUBLE_H

#include "psc.h"
#include "psc_particles_double.h"

#include "psc_fields_as_c.h"

typedef mparticles_double_t mparticles_t;
typedef struct psc_particles particles_t;
typedef particle_double_t particle_t;
typedef particle_double_real_t creal;
#define creal_abs(x) fabs(x)
#define creal_sqrt(x) sqrt(x)

#define psc_mparticles_get_cf    psc_mparticles_get_double
#define psc_mparticles_put_cf    psc_mparticles_put_double
#define psc_mparticles_get_patch psc_mparticles_get_patch_double
#define particles_get_one        particles_double_get_one
#define particle_qni_div_mni     particle_double_qni_div_mni
#define particle_qni_wni         particle_double_qni_wni
#define fint                     particle_double_real_fint

void psc_push_particles_double_1vb_push_yz(struct psc_push_particles *push,
					   mparticles_base_t *particles_base,
					   mfields_base_t *flds_base);

#endif
