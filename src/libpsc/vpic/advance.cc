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
#include "vpic_mfields.h"
#include "vpic_mparticles.h"
#include "vpic_push_particles.h"
#include "vpic_init.h"

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

void vpic_sort_run(struct vpic_mparticles *vmprts, int step)
{
  // Sort the particles for performance if desired.

  species_t *sp;

  LIST_FOR_EACH( sp, vmprts->species_list )
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
// vpic_push_particles

void vpic_push_particles::clear_accumulator_array()
{
  TIC ::clear_accumulator_array( accumulator_array ); TOC( clear_accumulators, 1 );
}

void vpic_push_particles::advance_p(vpic_mparticles *vmprts)
{
  species_t *sp;

  LIST_FOR_EACH( sp, vmprts->species_list )
    TIC ::advance_p( sp, accumulator_array, interpolator_array ); TOC( advance_p, 1 );
}

void vpic_push_particles::reduce_accumulator_array()
{
  TIC ::reduce_accumulator_array( accumulator_array ); TOC( reduce_accumulators, 1 );
}

void vpic_push_particles::boundary_p(vpic_mparticles *vmprts, vpic_mfields *vmflds)
{
  TIC
    for( int round=0; round<num_comm_round; round++ )
      ::boundary_p( simulation->particle_bc_list, vmprts->species_list,
		    vmflds, accumulator_array );
  TOC( boundary_p, num_comm_round );

  species_t *sp;
  LIST_FOR_EACH( sp, vmprts->species_list ) {
    if( sp->nm ) // && simulation->verbose )
      WARNING(( "Removing %i particles associated with unprocessed %s movers (increase num_comm_round)",
                sp->nm, sp->name ));
    // Drop the particles that have unprocessed movers due to a user defined
    // boundary condition. Particles of this type with unprocessed movers are
    // in the list of particles and move_p has set the voxel in the particle to
    // 8*voxel + face. This is an incorrect voxel index and in many cases can
    // in fact go out of bounds of the voxel indexing space. Removal is in
    // reverse order for back filling. Particle charge is accumulated to the
    // mesh before removing the particle.
    int nm = sp->nm;
    particle_mover_t * RESTRICT ALIGNED(16)  pm = sp->pm + sp->nm - 1;
    particle_t * RESTRICT ALIGNED(128) p0 = sp->p;
    for (; nm; nm--, pm--) {
      int i = pm->i; // particle index we are removing
      p0[i].i >>= 3; // shift particle voxel down
      // accumulate the particle's charge to the mesh
      accumulate_rhob( vmflds->f, p0+i, sp->g, sp->q );
      p0[i] = p0[sp->np-1]; // put the last particle into position i
      sp->np--; // decrement the number of particles
    }
    sp->nm = 0;
  }
}

void vpic_push_particles::unload_accumulator_array(vpic_mfields *vmflds)
{
  TIC ::unload_accumulator_array( vmflds, accumulator_array ); TOC( unload_accumulator, 1 );
}

void vpic_push_particles::load_interpolator_array(vpic_mfields *vmflds)
{
  TIC ::load_interpolator_array( interpolator_array, vmflds ); TOC( load_interpolator, 1 );
}

void vpic_push_particles::uncenter_p(vpic_mparticles *vmprts)
{
  species_t *sp;
  LIST_FOR_EACH( sp, vmprts->species_list )
    TIC ::uncenter_p( sp, interpolator_array ); TOC( uncenter_p, 1 );
}

// ======================================================================
// vpic_mfields

#define FAK kernel

void vpic_mfields::clear_jf()
{
  TIC FAK->clear_jf(this); TOC( clear_jf, 1 );
}

void vpic_mfields::synchronize_jf()
{
  TIC FAK->synchronize_jf(this); TOC( synchronize_jf, 1 );
}

void vpic_mfields::compute_div_b_err()
{
  TIC FAK->compute_div_b_err(this); TOC( compute_div_b_err, 1 );
}

void vpic_mfields::compute_div_e_err()
{
  TIC FAK->compute_div_e_err(this); TOC( compute_div_e_err, 1 );
}

double vpic_mfields::compute_rms_div_b_err()
{
  double err;
  TIC err = FAK->compute_rms_div_b_err(this); TOC( compute_rms_div_b_err, 1 );
  return err;
}

double vpic_mfields::compute_rms_div_e_err()
{
  double err;
  TIC err = FAK->compute_rms_div_e_err(this); TOC( compute_rms_div_e_err, 1 );
  return err;
}

void vpic_mfields::clean_div_b()
{
  TIC FAK->clean_div_b(this); TOC( clean_div_b, 1 );
}

void vpic_mfields::clean_div_e()
{
  TIC FAK->clean_div_e(this); TOC( clean_div_e, 1 );
}

void vpic_mfields::compute_curl_b()
{
  TIC FAK->compute_curl_b(this); TOC( compute_curl_b, 1 );
}

void vpic_mfields::clear_rhof()
{
  TIC FAK->clear_rhof(this); TOC( clear_rhof, 1 );
}

void vpic_mfields::accumulate_rho_p(vpic_mparticles *vmprts)
{
  species_t *sp;
  LIST_FOR_EACH( sp, vmprts->species_list )
    TIC ::accumulate_rho_p(this, sp); TOC( accumulate_rho_p, 1 );
}

void vpic_mfields::synchronize_rho()
{
  TIC FAK->synchronize_rho(this); TOC( synchronize_rho, 1 );
}

void vpic_mfields::compute_rhob()
{
  TIC FAK->compute_rhob(this); TOC( compute_rhob, 1 );
}

double vpic_mfields::synchronize_tang_e_norm_b()
{
  double err;
  TIC err = FAK->synchronize_tang_e_norm_b(this); TOC( synchronize_tang_e_norm_b, 1 );
  return err;
}

#undef FAK
#define FAK vmflds->kernel

// ======================================================================
// vpic_push_fields

void vpic_push_fields_advance_b(struct vpic_mfields *vmflds, double frac)
{
  TIC FAK->advance_b( vmflds, frac ); TOC( advance_b, 1 );
}

void vpic_push_fields_advance_e(struct vpic_mfields *vmflds, double frac)
{
  TIC FAK->advance_e( vmflds, frac ); TOC( advance_e, 1 );
}

// ======================================================================
// vpic_print_status

void vpic_print_status()
{
  update_profile( rank()==0 );
}

// ======================================================================
// vpic_moments

void vpic_moments_run(struct vpic_mfields_hydro *vmflds, struct vpic_mparticles *vmprts, int kind)
{
  // This relies on load_interpolator_array() having been called earlier
  clear_hydro_array(vmflds);
  species_t *sp;
  LIST_FOR_EACH(sp, vmprts->species_list) {
    if (sp->id == kind) {
      accumulate_hydro_p(vmflds, sp, simulation->interpolator_array);
      break;
    }
  }
  
  synchronize_hydro_array(vmflds);
}
