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

