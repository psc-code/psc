
#include "psc_push_particles_private.h"

#include "psc_push_fields_private.h"
#include "psc_bnd_fields.h"
#include "psc_bnd.h"
#include "push_particles.hxx"

#include <mrc_profile.h>

// ======================================================================
// forward to subclass

extern double *psc_balance_comp_time_by_patch;

void
psc_push_particles_prep(struct psc_push_particles *push,
			MparticlesBase& mprts_base, MfieldsBase& mflds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("push_particles_prep", 1., 0, 0);
  }  

  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);

  prof_start(pr);
  prof_restart(pr_time_step_no_comm);
  psc_stats_start(st_time_particle);

  if (ops->prep) {
    ops->prep(push, mprts_base, mflds_base);
  }

  psc_stats_stop(st_time_particle);
  prof_stop(pr_time_step_no_comm);
  prof_stop(pr);
}

// ======================================================================
// psc_push_particles class

struct mrc_class_psc_push_particles_ : mrc_class_psc_push_particles {
  mrc_class_psc_push_particles_() {
    name             = "psc_push_particles";
    size             = sizeof(struct psc_push_particles);
  }
} mrc_class_psc_push_particles;

