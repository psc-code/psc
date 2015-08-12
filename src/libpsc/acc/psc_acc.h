
#ifndef PSC_ACC_H
#define PSC_ACC_H

#include <psc_fields.h>
#include <psc_particles.h>

static const int psc_particles_acc_bs[3] = { 4, 4, 4 };

void acc_push_mflds_E_yz(struct psc_mfields *mflds);
void acc_push_mflds_H_yz(struct psc_mfields *mflds);

void acc_push_mflds_E_xyz(struct psc_mfields *mflds);
void acc_push_mflds_H_xyz(struct psc_mfields *mflds);

void acc_1vbec_push_mprts_yz(struct psc_mparticles *mprts,
			     struct psc_mfields *mflds);
void acc_1vbec_push_mprts_xyz(struct psc_mparticles *mprts,
			      struct psc_mfields *mflds);
#endif

