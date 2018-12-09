
// ======================================================================
// PushParticlesSimple

template<typename C, template<typename> class PushParticlesImpl>
struct PushParticlesSimple
{
  using Mparticles = typename C::Mparticles;
  using MfieldsState = typename C::MfieldsState;
  
  using checks_order = typename PushParticlesImpl<C>::checks_order;
  
  // ----------------------------------------------------------------------
  // push_mprts

  static void push_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    PushParticlesImpl<C>::push_mprts(mflds, mprts);
  }

  // ----------------------------------------------------------------------
  // stagger_mprts
  
  static void stagger_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    PushParticlesImpl<C>::stagger_mprts_patch(mflds, mprts);
  }
};

