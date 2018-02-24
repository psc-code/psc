
#include "psc_push_particles_private.h"

#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "psc_method.h"
#include "push_particles.hxx"
#include "vpic_iface.h"

// ======================================================================
// PushParticlesVpic

struct PushParticlesVpic : PushParticlesBase
{
  PushParticlesVpic()
  {
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);
  }

  void prep(struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base) override
  {
    // needs E, B
    mfields_vpic_t mf = mflds_base->get_as<mfields_vpic_t>(EX, HX + 6);
    FieldArray *vmflds = psc_mfields_vpic(mf.mflds())->vmflds_fields;
    
    Simulation_push_mprts_prep(sim_, vmflds);
    
    mf.put_as(mflds_base, 0, 0);
  }
  
  void push_mprts(struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base) override
  {
    // needs E, B (not really, because they're already in interpolator), rhob?
    mfields_vpic_t mf = mflds_base->get_as<mfields_vpic_t>(EX, HX + 6);
    FieldArray *vmflds = psc_mfields_vpic(mf.mflds())->vmflds_fields;
    mparticles_vpic_t mprts = mprts_base->get_as<mparticles_vpic_t>();
    
    Simulation_push_mprts(sim_, mprts->vmprts, vmflds);
    
    // update jf FIXME: rhob too, probably, depending on b.c.
    mf.put_as(mflds_base, JXI, JXI + 3);
    mprts.put_as(mprts_base);
  }

private:
  Simulation* sim_;
};

// ----------------------------------------------------------------------
// psc_push_particles_vpic_prep

static void
psc_push_particles_vpic_prep(struct psc_push_particles *push,
			     struct psc_mparticles *mprts_base,
			     struct psc_mfields *mflds_base)
{
  PscPushParticles<PushParticlesVpic> pushp(push);
  pushp->prep(mprts_base, mflds_base);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "vpic"

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticlesVpic>;

struct psc_push_particles_ops_vpic : psc_push_particles_ops {
  psc_push_particles_ops_vpic() {
    name                  = "vpic";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    push_mprts            = PushParticlesWrapper_t::push_mprts;
    prep                  = psc_push_particles_vpic_prep;
  }
} psc_push_particles_vpic_ops;

