
#ifndef PSC_PUSH_PARTICLES_1VB_H
#define PSC_PUSH_PARTICLES_1VB_H

#include <psc_push_particles_private.h>

void psc_push_particles_1vbec_single_by_block_push_mprts_yz(struct psc_push_particles *push,
							    struct psc_mparticles *mprts,
							    struct psc_mfields *mflds);
void psc_push_particles_1vbec_single_by_block_push_mprts_xyz(struct psc_push_particles *push,
							     struct psc_mparticles *mprts,
							     struct psc_mfields *mflds);

void psc_push_particles_1vbec_single_push_mprts_xyz(struct psc_push_particles *push,
						    struct psc_mparticles *mprts,
						    struct psc_mfields *mflds);
void psc_push_particles_1vbec_single_push_mprts_xyz_xz(struct psc_push_particles *push,
						       struct psc_mparticles *mprts,
						       struct psc_mfields *mflds);
void psc_push_particles_1vbec_single_push_mprts_yz(struct psc_push_particles *push,
						   struct psc_mparticles *mprts,
						   struct psc_mfields *mflds);
void psc_push_particles_1vbec_single_push_mprts_1(struct psc_push_particles *push,
						  struct psc_mparticles *mprts,
						  struct psc_mfields *mflds);

void psc_push_particles_1vbec_single_stagger_mprts_yz(struct psc_push_particles *push,
						      struct psc_mparticles *mprts,
						      struct psc_mfields *mflds);

void psc_push_particles_1vbec_double_push_mprts_yz(struct psc_push_particles *push,
						   struct psc_mparticles *mprts,
						   struct psc_mfields *mflds);
void psc_push_particles_1vbec_double_push_mprts_xyz(struct psc_push_particles *push,
						    struct psc_mparticles *mprts,
						    struct psc_mfields *mflds);
void psc_push_particles_1vbec_double_push_mprts_1(struct psc_push_particles *push,
						  struct psc_mparticles *mprts,
						  struct psc_mfields *mflds);
void psc_push_particles_1vbec_double_stagger_mprts_1(struct psc_push_particles *push,
						     struct psc_mparticles *mprts,
						     struct psc_mfields *mflds);


template<typename MF, typename dim>
struct push_particles_ops
{
  static void stagger_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			    struct psc_mfields *mflds_base);
};

#endif
