
#pragma once

#include "psc_push_particles_private.h"
#include "particles.hxx"
#include "fields3d.hxx"

// ======================================================================
// PscPushParticles

template<typename S>
struct PscPushParticles
{
  using sub_t = S;

  explicit PscPushParticles(psc_push_particles *pushp)
    : pushp_(pushp),
      sub_(mrc_to_subobj(pushp, sub_t))
  {}

  void operator()(mparticles_base_t mprts, mfields_base_t mflds)
  {
    psc_push_particles_run(pushp_, mprts.mprts(), mflds.mflds());
  }

  psc_push_particles *pushp() { return pushp_; }
  
  sub_t* operator->() { return sub_; }

  sub_t* sub() { return sub_; }

private:
  psc_push_particles *pushp_;
  sub_t *sub_;
};

// ======================================================================
// PushParticlesBase

class PushParticlesBase
{
public:
  virtual void prep(struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base)
  { assert(0); }
  
  virtual void push_mprts(struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void push_mprts_xyz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void push_mprts_xy(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void push_mprts_x(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void push_mprts_y(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void push_mprts_z(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void push_mprts_1(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { assert(0); }

  virtual void stagger_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { mprintf("WARNING: %s not implemented\n", __func__); }

  virtual void stagger_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { mprintf("WARNING: %s not implemented\n", __func__); }

  virtual void stagger_mprts_1(struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  { mprintf("WARNING: %s not implemented\n", __func__); }

};

using PscPushParticlesBase = PscPushParticles<PushParticlesBase>;

template<typename PushParticles_t>
class PushParticlesWrapper
{
public:
  const static size_t size = sizeof(PushParticles_t);
  
  static void setup(struct psc_push_particles *push)
  {
    PscPushParticles<PushParticles_t> pushp(push);
    new(pushp.sub()) PushParticles_t;
  }

  static void destroy(struct psc_push_particles *push)
  {
    PscPushParticles<PushParticles_t> pushp(push);
    pushp.sub()->~PushParticles_t();
  }

  static void push_mprts(struct psc_push_particles *push,
			 struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  {
    PscPushParticles<PushParticles_t> pushp(push);
    pushp->push_mprts(mprts, mflds_base);
  }

};

