/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 */

#include "vpic_iface.h"
#include "vpic_push_particles.h"
#include "field_array.h"
#include "hydro_array.h"

#include "vpic.h"

#include <cassert>

// ======================================================================

extern vpic_simulation *simulation;

int rank()
{
  return simulation->rank();
}

// ======================================================================
// vpic_sort

void vpic_sort_run(Particles *vmprts, int step)
{
  // Sort the particles for performance if desired.

  species_t *sp;

  LIST_FOR_EACH( sp, vmprts->sl_ )
    if( (sp->sort_interval>0) && ((step % sp->sort_interval)==0) ) {
      if( rank()==0 ) MESSAGE(( "Performance sorting \"%s\"", sp->name ));
      TIC sort_p( sp ); TOC( sort_p, 1 );
    } 
}

// ======================================================================
// vpic_push_fields

void vpic_push_fields_advance_b(FieldArray *vmflds, double frac)
{
  TIC vmflds->advance_b(frac); TOC(advance_b, 1);
}

void vpic_push_fields_advance_e(FieldArray *vmflds, double frac)
{
  TIC vmflds->advance_e(frac); TOC(advance_e, 1);
}

// ======================================================================
// vpic_moments

void vpic_moments_run(HydroArray *vmflds, Particles *vmprts, int kind)
{
  // This relies on load_interpolator_array() having been called earlier
  clear_hydro_array(vmflds);
  species_t *sp;
  LIST_FOR_EACH(sp, vmprts->sl_) {
    if (sp->id == kind) {
      accumulate_hydro_p(vmflds, sp, simulation->interpolator_array);
      break;
    }
  }
  
  synchronize_hydro_array(vmflds);
}
