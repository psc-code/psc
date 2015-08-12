
#include "psc_push_particles_private.h"

#include "psc_acc.h"
#include "psc_particles_acc.h"

#include <mrc_profile.h>

#include <string.h>

// ======================================================================
// psc_push_particles: subclass "1vb_acc"

// ----------------------------------------------------------------------
// psc_push_particles_1vbec_push_mprts_yz

static void
psc_push_particles_1vbec_push_mprts_yz(struct psc_push_particles *push,
				       struct psc_mparticles *mprts_base,	
				       struct psc_mfields *mflds_base)					
{									
  static int pr;
  if (!pr) {
    pr = prof_register("1vbec_yz", 1., 0, 0);
  }  

  struct psc_mparticles *mprts =
    psc_mparticles_get_as(mprts_base, "acc", 0);
  struct psc_mfields *mflds =
    psc_mfields_get_as(mflds_base, "acc", EX, EX + 6);

  prof_start(pr);
  acc_1vbec_push_mprts_yz(mprts, mflds);
  prof_stop(pr);

  psc_mparticles_put_as(mprts, mprts_base, 0);
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);	
}

// ----------------------------------------------------------------------
// psc_push_particles_1vbec_push_mprts_xyz

static void
psc_push_particles_1vbec_push_mprts_xyz(struct psc_push_particles *push,
					struct psc_mparticles *mprts_base,	
					struct psc_mfields *mflds_base)					
{									
  static int pr;
  if (!pr) {
    pr = prof_register("1vbec_xyz", 1., 0, 0);
  }  
  static int pr_2;
  if (!pr_2) {
    pr_2 = prof_register("1vbec_xyz 2", 1., 0, 0);
  }  

  struct psc_mparticles *mprts =
    psc_mparticles_get_as(mprts_base, "acc", 0);
  struct psc_mfields *mflds =
    psc_mfields_get_as(mflds_base, "acc", EX, EX + 6);

  prof_start(pr);
  acc_1vbec_push_mprts_xyz(mprts, mflds);
  prof_stop(pr);

  psc_mparticles_put_as(mprts, mprts_base, 0);
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);	
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vbec_acc"
									
struct psc_push_particles_ops						
psc_push_particles_1vbec_acc_ops = {	
  .name                  = "1vbec_acc",
  .push_mprts_yz         = psc_push_particles_1vbec_push_mprts_yz,
  .push_mprts_xyz        = psc_push_particles_1vbec_push_mprts_xyz,
};
