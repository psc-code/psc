
#pragma once

#include "push_particles.hxx"
#include "vpic_iface.h"

// ======================================================================
// PushParticlesVpic

struct PushParticlesVpic : PushParticlesBase
{
  using Mparticles = MparticlesVpic;
  using MfieldsState = MfieldsState;
  
  void push_mprts(Mparticles& mprts, MfieldsState& mflds, MfieldsInterpolator& interpolator,
		  MfieldsAccumulator& accumulator, ParticleBcList& particle_bc_list,
		  int num_comm_round)
  {
    // For this to work, interpolator needs to have been set from mflds E/B before,
    // ie., we're not using mflds for E and B here at all.
    
    // At this point, fields are at E_0 and B_0 and the particle positions
    // are at r_0 and u_{-1/2}.  Further the mover lists for the particles should
    // empty and all particles should be inside the local computational domain.
    // Advance the particle lists.
    if (!mprts.empty()) {
      TIC AccumulatorOps::clear(accumulator); TOC(clear_accumulators, 1);
      ParticlesOps::advance_p(mprts, accumulator, interpolator);
    }

    // Because the partial position push when injecting aged particles might
    // place those particles onto the guard list (boundary interaction) and
    // because advance_p requires an empty guard list, particle injection must
    // be done after advance_p and before guard list processing. Note:
    // user_particle_injection should be a stub if sl_ is empty.
#if 0
    if (emitter_list) {
      TIC ::apply_emitter_list(emitter_list); TOC(emission_model, 1);
    }
    TIC user_particle_injection(); TOC(user_particle_injection, 1);
#endif
    
    // This should be after the emission and injection to allow for the
    // possibility of thread parallelizing these operations
    if (!mprts.empty()) {
      TIC AccumulatorOps::reduce(accumulator); TOC(reduce_accumulators, 1);
    }
    
    // At this point, most particle positions are at r_1 and u_{1/2}. Particles
    // that had boundary interactions are now on the guard list. Process the
    // guard lists. Particles that absorbed are added to rhob (using a corrected
    // local accumulation).
    TIC
      for(int round = 0; round < num_comm_round; round++) {
	ParticlesOps::boundary_p(particle_bc_list, mprts, mflds, accumulator);
      } TOC(boundary_p, num_comm_round);
    
    // Drop the particles that have unprocessed movers at this point
    ParticlesOps::drop_p(mprts, mflds);

    // At this point, all particle positions are at r_1 and u_{1/2}, the
    // guard lists are empty and the accumulators on each processor are current.
    // Convert the accumulators into currents.
    TIC AccumulateOps::clear_jf(mflds); TOC(clear_jf, 1);
    if (!mprts.empty()) {
      TIC AccumulatorOps::unload(accumulator, mflds); TOC(unload_accumulator, 1);
    }
    TIC AccumulateOps::synchronize_jf(mflds); TOC(synchronize_jf, 1);

    // At this point, the particle currents are known at jf_{1/2}.
    // Let the user add their own current contributions. It is the users
    // responsibility to insure injected currents are consistent across domains.
    // It is also the users responsibility to update rhob according to
    // rhob_1 = rhob_0 + div juser_{1/2} (corrected local accumulation) if
    // the user wants electric field divergence cleaning to work.
#if 0
    TIC user_current_injection(); TOC(user_current_injection, 1);
#endif
  }

  void prep(Mparticles& mprts, MfieldsState& mflds, MfieldsInterpolator& interpolator)
  {
    // At end of step:
    // Fields are updated ... load the interpolator for next time step and
    // particle diagnostics in user_diagnostics if there are any particle
    // species to worry about
    
    if (!mprts.empty()) {
      TIC InterpolatorOps::load(interpolator, mflds); TOC(load_interpolator, 1);
    }
  }
  
  void prep(MparticlesBase& mprts_base, MfieldsStateBase& mflds_base) override
  {
    assert(0);
#if 0
    // needs E, B
    auto mprts = mprts_base.get_as<Mparticles>();
    auto& mflds = mflds_base.get_as<Mfields>(EX, HX + 6);
    prep(mprts, mflds);
    mflds_base.put_as(mflds, 0, 0);
    mprts_base.put_as(mprts, MP_DONT_COPY);
#endif
  }
  
  void push_mprts(MparticlesBase& mprts_base, MfieldsStateBase& mflds_base) override
  {
    assert(0);
#if 0
    // needs E, B (not really, because they're already in interpolator), rhob?
    auto& mflds = mflds_base.get_as<Mfields>(EX, HX + 6);
    auto& mprts = mprts_base.get_as<Mparticles>();

    push_mprts(mprts, mflds);
    
    // update jf FIXME: rhob too, probably, depending on b.c.
    mflds_base.put_as(mflds, JXI, JXI + 3);
    mprts_base.put_as(mprts);
#endif
  }
};

