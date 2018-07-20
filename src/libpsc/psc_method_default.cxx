
#include "psc_method_private.h"

#include <psc_diag.h>
#include <psc_output_fields_collection.h>
#include <psc_output_particles.h>
#include <psc_particles_single.h>
#include <psc_particles_double.h>
#include <particles.hxx>
#include <push_particles.hxx>
#include <setup_particles.hxx>
#include <output_particles.hxx>

#include <stdlib.h>

// ======================================================================
// psc_method "default"

struct psc_method_ops_default : psc_method_ops {
  psc_method_ops_default() {
    name                          = "default";
  }
} psc_method_ops_default;
