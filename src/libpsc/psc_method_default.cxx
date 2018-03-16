
#include "psc_method_private.h"

#include <psc_push_particles.h>
#include <psc_diag.h>
#include <psc_output_fields_collection.h>
#include <psc_output_particles.h>
#include <psc_balance.h>
#include <psc_particles_single.h>
#include <psc_particles_double.h>
#include <particles.hxx>
#include <push_particles.hxx>
#include <setup_particles.hxx>

#include <stdlib.h>

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
				   uint *n_prts_by_patch)
{
  SetupParticles<MparticlesDouble>::setup_partition(psc, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_method_default_set_ic_particles

static void
psc_method_default_set_ic_particles(struct psc_method *method, struct psc *psc,
				    uint *n_prts_by_patch)
{
  PscMparticlesBase mprts(psc->particles);
  mprts->reserve_all(n_prts_by_patch);
  SetupParticles<MparticlesDouble>::setup_particles(psc, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_method_default_set_ic_fields

static void
psc_method_default_set_ic_fields(struct psc_method *method, struct psc *psc)
{
  // set fields E^{n+1/2}, B^{n+1/2}
  psc_set_ic_fields(psc);
}

// ----------------------------------------------------------------------
// psc_method_default_initialize

static void
psc_method_default_initialize(struct psc_method *method, struct psc *psc)
{
  auto pushp = PscPushParticlesBase{psc->push_particles};
  pushp.stagger(PscMparticlesBase{psc->particles}, PscMfieldsBase{psc->flds});

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

struct psc_method_ops_default : psc_method_ops {
  psc_method_ops_default() {
    name                          = "default";
    do_setup                      = psc_method_default_do_setup;
    setup_partition               = psc_method_default_setup_partition;
    set_ic_particles              = psc_method_default_set_ic_particles;
    set_ic_fields                 = psc_method_default_set_ic_fields;
    initialize                    = psc_method_default_initialize;
    output                        = psc_method_default_output;
  }
} psc_method_ops_default;
