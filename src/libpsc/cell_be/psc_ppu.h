#ifndef PSC_CBE_H
#define PSC_CBE_H

// Order of includes is important!
// There are preprocessor sections in 
// psc_cbe_common.h which depend on 
// which include files preceded it. 
// I have a feeling this is very bad form.
#include "psc.h"
#include "psc_particles_as_cbe.h"
#include "psc_fields_as_c.h"

#include "psc_cbe_common.h"

#ifdef CELLEMU
#include "libspe2_c.h"
#else
#include <libspe2.h>
#endif 

#ifdef CELLEMU
#define NR_SPE (1)
#else
#define NR_SPE (8)
#endif 


// Now that Kai is transitioning to the 'patch' multidomain structure, 
// the entire block structure I wrote is not needed and has been removed.
// This is a happy day, kiddies, cause that stuff was ugly.
//
// I also means we're going to drop the specialized sorting (which doesn't 
// actually work), and seriously cuts down on the amount of specialized ppu
// code we need to have. 

void psc_push_particles_cbe_push_xy(struct psc_push_particles *push, struct psc_mparticles *particles_base, struct psc_mfields *flds_base);


// Spe handeling functions from psc_cbe.c
void psc_init_spes(void);
void psc_kill_spes(void);
void update_spes_status(void);
void cell_run_patch(int p, struct psc_fields *pf, particles_t *pp, int job);
void wait_all_spe(void);

void cbe_create(void);
void cbe_destroy(void);

#endif
