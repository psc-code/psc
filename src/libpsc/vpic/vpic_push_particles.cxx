
#include "vpic_push_particles.h"

#include "field_array.h"

// ======================================================================
// vpic_push_particles implementation

// ----------------------------------------------------------------------
// ctor

vpic_push_particles::vpic_push_particles(Simulation *sim)
  : sim_(sim)
{
  interpolator = sim->interpolator_;
  accumulator = sim->accumulator_;
  num_comm_round = sim->num_comm_round_;
}

// ======================================================================
// wrappers

// ----------------------------------------------------------------------
// vpic_push_particles_new_from_Simulation

struct vpic_push_particles *
vpic_push_particles_new_from_Simulation(Simulation *sim)
{
  return new vpic_push_particles(sim);
}

// ----------------------------------------------------------------------
// vpic_push_particles forwards

void vpic_push_particles_push_mprts(struct vpic_push_particles *vpushp,
				    Particles *vmprts,
				    FieldArray *vmflds)
{
  vpushp->push_mprts(vmprts, vmflds);
}

void vpic_push_particles_prep(struct vpic_push_particles *vpushp,
			      Particles *vmprts, FieldArray *vmflds)
{
  vpushp->prep(vmprts, vmflds);
}

void vpic_push_particles_stagger_mprts(struct vpic_push_particles *vpushp,
				       Particles *vmprts, FieldArray *vmflds)
{
  vpushp->stagger_mprts(vmprts, vmflds);
}

// ----------------------------------------------------------------------

void vpic_push_particles::stagger_mprts(Particles *vmprts, FieldArray *vmflds)
{
  if (!vmprts->empty()) {
    sim_->load_interpolator_array(interpolator, vmflds);

    for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
      TIC ::uncenter_p(&*sp, interpolator); TOC(uncenter_p, 1);
    }
  }
}

void vpic_push_particles::push_mprts(Particles *vmprts, FieldArray *vmflds)
{
  // At this point, fields are at E_0 and B_0 and the particle positions
  // are at r_0 and u_{-1/2}.  Further the mover lists for the particles should
  // empty and all particles should be inside the local computational domain.
  // Advance the particle lists.

  if (!vmprts->empty()) {
    sim_->clear_accumulator_array(accumulator);
    sim_->advance_p(*vmprts, *accumulator, *interpolator);
  }

  // Because the partial position push when injecting aged particles might
  // place those particles onto the guard list (boundary interaction) and
  // because advance_p requires an empty guard list, particle injection must
  // be done after advance_p and before guard list processing. Note:
  // user_particle_injection should be a stub if sl_ is empty.

  sim_->emitter();

  if (!vmprts->empty()) {
    // This should be after the emission and injection to allow for the
    // possibility of thread parallelizing these operations

    sim_->reduce_accumulator_array(accumulator);
  }
  // At this point, most particle positions are at r_1 and u_{1/2}. Particles
  // that had boundary interactions are now on the guard list. Process the
  // guard lists. Particles that absorbed are added to rhob (using a corrected
  // local accumulation).

  TIC
    for( int round=0; round<num_comm_round; round++ )
      ::boundary_p( sim_->simulation_->particle_bc_list, vmprts->head(),
		    vmflds, accumulator );
  TOC( boundary_p, num_comm_round );

  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
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

#if 0
  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    for (int n = 0; n < sp->np; n++) {
      particle_t *p = &sp->p[n];
      int i = p->i;
      int im[3] = { sp->g->nx + 2, sp->g->ny + 2, sp->g->nz + 2 };
      int i3[3];
      i3[2] = i / (im[0] * im[1]); i -= i3[2] * (im[0] * im[1]);
      i3[1] = i / im[0]; i-= i3[1] * im[0];
      i3[0] = i;
      if (!(i3[2] >= 1 && i3[2] <= sp->g->nz)) {
	mprintf("i3 %d %d %d\n", i3[0], i3[1], i3[2]);
	assert(0);
      }
    }
  }
#endif
  
  // At this point, all particle positions are at r_1 and u_{1/2}, the
  // guard lists are empty and the accumulators on each processor are current.
  // Convert the accumulators into currents.

  vmflds->clear_jf();
  if (!vmprts->empty()) {
    sim_->unload_accumulator_array(vmflds, accumulator);
  }
  vmflds->synchronize_jf();

  // At this point, the particle currents are known at jf_{1/2}.
  // Let the user add their own current contributions. It is the users
  // responsibility to insure injected currents are consistent across domains.
  // It is also the users responsibility to update rhob according to
  // rhob_1 = rhob_0 + div juser_{1/2} (corrected local accumulation) if
  // the user wants electric field divergence cleaning to work.

  sim_->current_injection();
}

void vpic_push_particles::prep(Particles *vmprts, FieldArray *vmflds)
{
  if (!vmprts->empty()) {
    sim_->load_interpolator_array(interpolator, vmflds);
  }
}

