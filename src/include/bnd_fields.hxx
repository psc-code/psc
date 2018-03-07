
#pragma once

#include "psc_bnd_fields_private.h"

// ======================================================================
// BndFieldsBase

struct BndFieldsBase
{
  virtual void fill_ghosts(PscMfieldsBase mflds, int mb, int me) = 0;
  virtual void add_ghosts(PscMfieldsBase mflds, int mb, int me) = 0;
  virtual void reset() = 0;
};

// ======================================================================
// PscBndFields

template<typename S>
struct PscBndFields
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, BndFieldsBase*>::value,
  		"sub classes used in PscBndFieldsParticles must derive from BndFieldsBase");
  
  explicit PscBndFields(psc_bnd_fields *bndf)
    : bndf_(bndf)
  {}

  void fill_ghosts_E(PscMfieldsBase mflds)
  {
    psc_bnd_fields_fill_ghosts_E(bndf_, mflds.mflds());
  }

  void fill_ghosts_H(PscMfieldsBase mflds)
  {
    psc_bnd_fields_fill_ghosts_H(bndf_, mflds.mflds());
  }

  void add_ghosts_J(PscMfieldsBase mflds)
  {
    psc_bnd_fields_add_ghosts_J(bndf_, mflds.mflds());
  }

  sub_t* sub() { return mrc_to_subobj(bndf_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_bnd_fields* bndf_;
};

using PscBndFieldsBase = PscBndFields<BndFieldsBase>;

