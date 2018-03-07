
#pragma once

#include "psc_bnd_fields_private.h"

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
    struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(bndf_);
    if (ops->fill_ghosts_E) {
      ops->fill_ghosts_E(bndf_, mflds.mflds());
    }
  }

  void fill_ghosts_H(PscMfieldsBase mflds)
  {
    struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(bndf_);
    if (ops->fill_ghosts_H) {
      ops->fill_ghosts_H(bndf_, mflds.mflds());
    }
  }

  void add_ghosts_J(PscMfieldsBase mflds)
  {
    struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(bndf_);
    if (ops->add_ghosts_J) {
      ops->add_ghosts_J(bndf_, mflds.mflds());
    }
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

  static void fill_ghosts_E(psc_bnd_fields *_bndf, psc_mfields *mflds_base)
  {
    PscBndFields<BndFields> bndf(_bndf);
    bndf->fill_ghosts_E(mflds_base);
  }
  
  static void fill_ghosts_H(psc_bnd_fields *_bndf, psc_mfields *mflds_base)
  {
    PscBndFields<BndFields> bndf(_bndf);
    bndf->fill_ghosts_H(mflds_base);
  }
  
  static void add_ghosts_J(psc_bnd_fields *_bndf, psc_mfields *mflds_base)
  {
    PscBndFields<BndFields> bndf(_bndf);
    bndf->add_ghosts_J(mflds_base);
  }
  
};

