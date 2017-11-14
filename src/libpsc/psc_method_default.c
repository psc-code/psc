
#include "psc_method_private.h"

#include <psc_push_particles.h>
#include <psc_diag.h>
#include <psc_output_fields_collection.h>
#include <psc_output_particles.h>

// ======================================================================
// psc_method "default"

// ----------------------------------------------------------------------
// psc_method_default_do_setup

static void
psc_method_default_do_setup(struct psc_method *method, struct psc *psc)
{
  psc_setup_coeff(psc);
  psc_setup_domain(psc);
}

// ----------------------------------------------------------------------
// psc_method_default_setup_partition

static void
psc_method_default_setup_partition(struct psc_method *method, struct psc *psc,
				   int *n_prts_by_patch, int *particle_label_offset)
{
  psc_setup_partition(psc, n_prts_by_patch, particle_label_offset);
}

// ----------------------------------------------------------------------
// psc_method_default_setup_partition_and_particles

static void
psc_method_default_setup_partition_and_particles(struct psc_method *method, struct psc *psc)
{
  // alloc / initialize particles
  int particle_label_offset;
  int *nr_particles_by_patch = calloc(psc->nr_patches, sizeof(*nr_particles_by_patch));
  psc_method_setup_partition(psc->method, psc, nr_particles_by_patch, &particle_label_offset);
  psc_balance_initial(psc->balance, psc, &nr_particles_by_patch);

  psc->particles = 
    psc_mparticles_create(mrc_domain_comm(psc->mrc_domain));
  psc_mparticles_set_type(psc->particles, psc->prm.particles_base);
  psc_mparticles_set_name(psc->particles, "mparticles");
  int nr_patches;
  mrc_domain_get_patches(psc->mrc_domain, &nr_patches);
  psc_mparticles_set_param_int(psc->particles, "nr_patches", nr_patches);
  if (psc->prm.particles_base_flags == 0) {
    psc->prm.particles_base_flags = psc_push_particles_get_mp_flags(ppsc->push_particles);
  }
  psc_mparticles_set_param_int(psc->particles, "flags", psc->prm.particles_base_flags);
  psc_mparticles_setup(psc->particles);
  psc_mparticles_reserve_all(psc->particles, nr_particles_by_patch);

  psc_setup_particles(psc, nr_particles_by_patch, particle_label_offset);
  free(nr_particles_by_patch);
}

// ----------------------------------------------------------------------
// psc_method_default_setup_fields

static void
psc_method_default_setup_fields(struct psc_method *method, struct psc *psc)
{
  // set fields E^{n+1/2}, B^{n+1/2}
  psc_setup_fields(psc);
}

// ----------------------------------------------------------------------
// psc_method_default_initialize

static void
psc_method_default_initialize(struct psc_method *method, struct psc *psc)
{
  psc_push_particles_stagger(psc->push_particles, psc->particles, psc->flds);

  // initial output / stats
  psc_output(psc);
  psc_stats_log(psc);
  psc_print_profiling(psc);
}

// ----------------------------------------------------------------------
// psc_method_default_output

void
psc_method_default_output(struct psc_method *method, struct psc *psc)
{
  psc_diag_run(psc->diag, psc);
  psc_output_fields_collection_run(psc->output_fields_collection, psc->flds, psc->particles);
  psc_output_particles_run(psc->output_particles, psc->particles);
}

// ----------------------------------------------------------------------
// psc_method "default"

struct psc_method_ops psc_method_ops_default = {
  .name                          = "default",
  .do_setup                      = psc_method_default_do_setup,
  .setup_fields                  = psc_method_default_setup_fields,
  .setup_partition               = psc_method_default_setup_partition,
  .setup_partition_and_particles = psc_method_default_setup_partition_and_particles,
  .initialize                    = psc_method_default_initialize,
  .output                        = psc_method_default_output,
};
