
#ifndef PSC_PUSH_PARTICLES_SINGLE_H
#define PSC_PUSH_PARTICLES_SINGLE_H

#include "psc.h"
#include "psc_particles_single.h"

#include "psc_fields_as_c.h"

typedef mparticles_single_t mparticles_t;
typedef particles_single_t particles_t;
typedef particle_single_t particle_t;
typedef particle_single_real_t creal;
#define creal_abs(x) fabsf(x)
#define creal_sqrt(x) sqrtf(x)

#define psc_mparticles_get_cf    psc_mparticles_get_single
#define psc_mparticles_put_cf    psc_mparticles_put_single
#define psc_mparticles_get_patch psc_mparticles_get_patch_single
#define particles_get_one        particles_single_get_one
#define particle_qni_div_mni     particle_single_qni_div_mni
#define particle_qni_wni         particle_single_qni_wni
#define fint                     particle_single_real_fint

void psc_push_particles_single_1vb_push_yz(struct psc_push_particles *push,
					   mparticles_base_t *particles_base,
					   mfields_base_t *flds_base);

#endif
