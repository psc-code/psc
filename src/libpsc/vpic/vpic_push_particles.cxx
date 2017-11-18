
#include "vpic_push_particles.h"

#include "vpic_mparticles.h"
#include "field_array.h"

extern vpic_simulation *simulation;

// ======================================================================
// vpic_push_particles

// ----------------------------------------------------------------------
// vpic_push_particles_create

struct vpic_push_particles *
vpic_push_particles_create()
{
  return new vpic_push_particles;
}

// ----------------------------------------------------------------------
// vpic_push_particles_ctor_from_simulation

void
vpic_push_particles_ctor_from_simulation(struct vpic_push_particles *vpushp,
					 struct Simulation *sim)
{
  vpushp->interpolator_array = static_cast<InterpolatorArray*>(sim->interpolator_array_);
  vpushp->accumulator_array = static_cast<AccumulatorArray*>(sim->accumulator_array_);
  vpushp->num_comm_round = simulation->num_comm_round;
}

// ----------------------------------------------------------------------
// vpic_push_particles_push_mprts

void vpic_push_particles_push_mprts(struct vpic_push_particles *vpushp,
				    struct vpic_mparticles *vmprts,
				    struct FieldArray *vmflds)
{
  // FIXME, this is kinda too much stuff all in here,
  // so it should be split up, but it'll do for now

  // At this point, fields are at E_0 and B_0 and the particle positions
  // are at r_0 and u_{-1/2}.  Further the mover lists for the particles should
  // empty and all particles should be inside the local computational domain.
  // Advance the particle lists.

  if (vmprts->sl_) {
    vpushp->clear_accumulator_array();
    vpushp->advance_p(vmprts);
  }

  // Because the partial position push when injecting aged particles might
  // place those particles onto the guard list (boundary interaction) and
  // because advance_p requires an empty guard list, particle injection must
  // be done after advance_p and before guard list processing. Note:
  // user_particle_injection should be a stub if sl_ is empty.

  vpic_emitter();

  if (vmprts->sl_) {
    // This should be after the emission and injection to allow for the
    // possibility of thread parallelizing these operations

    vpushp->reduce_accumulator_array();
  }
  // At this point, most particle positions are at r_1 and u_{1/2}. Particles
  // that had boundary interactions are now on the guard list. Process the
  // guard lists. Particles that absorbed are added to rhob (using a corrected
  // local accumulation).

  vpushp->boundary_p(vmprts, vmflds);

  // At this point, all particle positions are at r_1 and u_{1/2}, the
  // guard lists are empty and the accumulators on each processor are current.
  // Convert the accumulators into currents.

  vmflds->clear_jf();
  if (vmprts->sl_) {
    vpushp->unload_accumulator_array(vmflds);
  }
  vmflds->synchronize_jf();

  // At this point, the particle currents are known at jf_{1/2}.
  // Let the user add their own current contributions. It is the users
  // responsibility to insure injected currents are consistent across domains.
  // It is also the users responsibility to update rhob according to
  // rhob_1 = rhob_0 + div juser_{1/2} (corrected local accumulation) if
  // the user wants electric field divergence cleaning to work.

  vpic_current_injection();
}

// ----------------------------------------------------------------------
// vpic_push_particles_prep

void vpic_push_particles_prep(struct vpic_push_particles *vpushp,
			      struct vpic_mparticles *vmprts, struct FieldArray *vmflds)
{
  if (vmprts->sl_) {
    vpushp->load_interpolator_array(vmflds);
  }
}

// ----------------------------------------------------------------------
// vpic_push_particles_stagger_mprts

void vpic_push_particles_stagger_mprts(struct vpic_push_particles *vpushp,
				       struct vpic_mparticles *vmprts,
				       struct FieldArray *vmflds)
{
  if (vmprts->sl_) {
    vpushp->load_interpolator_array(vmflds);
  }
  vpushp->uncenter_p(vmprts);
}
