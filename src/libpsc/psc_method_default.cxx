
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
#include <output_particles.hxx>

#include <stdlib.h>

// ======================================================================
// psc_method "default"

// ----------------------------------------------------------------------
// psc_method_default_output

void
psc_method_default_output(struct psc_method *method, struct psc *psc, MparticlesBase& mprts)
{
  psc_diag_run(psc->diag, psc, mprts);
  psc_output_fields_collection_run(psc->output_fields_collection, psc->flds, mprts);
  PscOutputParticlesBase{psc->output_particles}.run(mprts);
}

// ----------------------------------------------------------------------
// psc_method "default"

struct psc_method_ops_default : psc_method_ops {
  psc_method_ops_default() {
    name                          = "default";
    output                        = psc_method_default_output;
  }
} psc_method_ops_default;
