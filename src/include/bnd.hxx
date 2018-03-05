
#pragma once

#include "psc_bnd_private.h"
#include "fields3d.hxx"

#include <mrc_profile.h>

// ======================================================================
// BndBase

struct BndBase
{
  virtual void fill_ghosts() = 0;
};

// ======================================================================
// PscBnd

template<typename S>
struct PscBnd
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, BndBase*>::value,
  		"sub classes used in PscBndParticles must derive from BndBase");
  
  explicit PscBnd(psc_bnd *bnd)
    : bnd_(bnd)
  {}

  void fill_ghosts(PscMfieldsBase mflds, int mb, int me)
  {
    psc_bnd_fill_ghosts(bnd_, mflds.mflds(), mb, me);
  }

  void add_ghosts(PscMfieldsBase mflds, int mb, int me)
  {
    psc_bnd_add_ghosts(bnd_, mflds.mflds(), mb, me);
  }

  sub_t* sub() { return mrc_to_subobj(bnd_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_bnd* bnd_;
};

using PscBndBase = PscBnd<BndBase>;

