
#pragma once

#include "psc_bnd_fields_private.h"
#include "fields3d.hxx"

// ======================================================================
// BndFieldsBase

struct BndFieldsBase
{
  virtual void fill_ghosts_E(PscMfieldsBase mflds_base) = 0;
  virtual void fill_ghosts_H(PscMfieldsBase mflds_base) = 0;
  virtual void add_ghosts_J(PscMfieldsBase mflds_base) = 0;
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
    sub()->fill_ghosts_E(mflds);
  }

  void fill_ghosts_H(PscMfieldsBase mflds)
  {
    sub()->fill_ghosts_H(mflds);
  }

  void add_ghosts_J(PscMfieldsBase mflds)
  {
    sub()->add_ghosts_J(mflds);
  }

  sub_t* sub() { return mrc_to_subobj(bndf_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_bnd_fields* bndf_;
};

using PscBndFieldsBase = PscBndFields<BndFieldsBase>;

// ======================================================================
// PscBndWrapper

template<typename BndFields>
class PscBndFieldsWrapper
{
public:
  const static size_t size = sizeof(BndFields);
  
  static void setup(psc_bnd_fields* _bndf)
  {
    PscBndFields<BndFields> bndf(_bndf);
    new(bndf.sub()) BndFields();
  }

  static void destroy(psc_bnd_fields* _bndf)
  {
    PscBndFields<BndFields> bndf(_bndf);
    bndf->~BndFields();
  }
};

