
#include "psc_bnd_fields_private.h"
#include "bnd_fields.hxx"

// FIXME, this is useful at most for testing and maybe should go away

struct BndFieldsNone : BndFieldsBase
{
  void fill_ghosts_E(PscMfieldsBase mflds_base) override {};
  void fill_ghosts_H(PscMfieldsBase mflds_base) override {};
  void add_ghosts_J(PscMfieldsBase mflds_base) override {};
};

// ======================================================================
// psc_bnd_fields: subclass "none"

struct psc_bnd_fields_ops_none : psc_bnd_fields_ops {
  using PscBndFields_t = PscBndFieldsWrapper<BndFieldsNone>;
  psc_bnd_fields_ops_none() {
    name                  = "none";
    size                  = PscBndFields_t::size;
    setup                 = PscBndFields_t::setup;
    destroy               = PscBndFields_t::destroy;
  }
} psc_bnd_fields_none_ops;
