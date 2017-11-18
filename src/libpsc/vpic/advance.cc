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
// vpic_collision

void vpic_collision_run()
{
  // Note: Particles should not have moved since the last performance sort
  // when calling collision operators.
  // FIXME: Technically, this placement of the collision operators only
  // yields a first order accurate Trotter factorization (not a second
  // order accurate factorization).

  if( simulation->collision_op_list ) {
    // FIXME: originally, vpic_clear_accumulator_array() was called before this.
    // It's now called later, though. I'm not sure why that would be necessary here,
    // but it needs to be checked.
    // The assert() below doesn't unfortunately catch all cases where this might go wrong
    // (ie., it's missing the user_particle_collisions())

    assert(0);
    TIC apply_collision_op_list( simulation->collision_op_list ); TOC( collision_model, 1 );
  }
  TIC simulation->user_particle_collisions(); TOC( user_particle_collisions, 1 );
}

// ======================================================================
// vpic_emitter

void vpic_emitter()
{
  if( simulation->emitter_list )
    TIC apply_emitter_list( simulation->emitter_list ); TOC( emission_model, 1 );
  TIC simulation->user_particle_injection(); TOC( user_particle_injection, 1 );
}

// ======================================================================
// vpic_current_injection

void vpic_current_injection()
{
  TIC simulation->user_current_injection(); TOC( user_current_injection, 1 );
}

// ======================================================================
// vpic_field_injection

void vpic_field_injection()
{
  // Let the user add their own contributions to the electric field. It is the
  // users responsibility to insure injected electric fields are consistent
  // across domains.

  TIC simulation->user_field_injection(); TOC( user_field_injection, 1 );
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
