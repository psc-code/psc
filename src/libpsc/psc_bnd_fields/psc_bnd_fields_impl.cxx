
#include "psc_bnd_fields_private.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"

#include "psc_bnd_fields_impl.hxx"

// ======================================================================
// psc_bnd_fields: subclass "c"


struct psc_bnd_fields_ops_c : psc_bnd_fields_ops {
  using BndFields_t = BndFields_<MfieldsC>;
  using PscBndFields_t = PscBndFieldsWrapper<BndFields_t>;
  psc_bnd_fields_ops_c() {
    name                  = "c";
    size                  = PscBndFields_t::size;
    setup                 = PscBndFields_t::setup;
    destroy               = PscBndFields_t::destroy;
  }
} psc_bnd_fields_c_ops;

// ======================================================================
// psc_bnd_fields: subclass "single"

struct psc_bnd_fields_ops_single : psc_bnd_fields_ops {
  using BndFields_t = BndFields_<MfieldsSingle>;
  using PscBndFields_t = PscBndFieldsWrapper<BndFields_t>;
  psc_bnd_fields_ops_single() {
    name                  = "single";
    size                  = PscBndFields_t::size;
    setup                 = PscBndFields_t::setup;
    destroy               = PscBndFields_t::destroy;
  }
} psc_bnd_fields_single_ops;

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
