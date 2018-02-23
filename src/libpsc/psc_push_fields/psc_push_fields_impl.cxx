
#include "psc_push_fields_private.h"

#include "psc_fields_single.h"
#include "psc_fields_c.h"

#include "psc_push_fields_impl.hxx"

// ----------------------------------------------------------------------
// psc_push_fields_sub_push_mflds_E
//

template<typename mfields_t>
static void psc_push_fields_sub_push_mflds_E(struct psc_push_fields *push,
					     struct psc_mfields *mflds_base,
					     double dt_fac)
{
  PushFields<mfields_t> pushf;
  pushf.push_E(push, mflds_base, dt_fac);
}

// ----------------------------------------------------------------------
// psc_push_fields_sub_push_mflds_H

template<typename mfields_t>
static void psc_push_fields_sub_push_mflds_H(struct psc_push_fields *push,
					     struct psc_mfields *mflds_base,
					     double dt_fac)
{
  PushFields<mfields_t> pushf;
  pushf.push_H(push, mflds_base, dt_fac);
}

// ======================================================================
// psc_push_fields: subclass "single"

struct psc_push_fields_ops_single : psc_push_fields_ops {
  psc_push_fields_ops_single() {
    name                  = "single";
    push_mflds_E          = psc_push_fields_sub_push_mflds_E<mfields_single_t>;
    push_mflds_H          = psc_push_fields_sub_push_mflds_H<mfields_single_t>;
  }
} psc_push_fields_single_ops;

// ======================================================================
// psc_push_fields: subclass "c"

struct psc_push_fields_ops_c : psc_push_fields_ops {
  psc_push_fields_ops_c() {
    name                  = "c";
    push_mflds_E          = psc_push_fields_sub_push_mflds_E<mfields_c_t>;
    push_mflds_H          = psc_push_fields_sub_push_mflds_H<mfields_c_t>;
  }
} psc_push_fields_c_ops;
