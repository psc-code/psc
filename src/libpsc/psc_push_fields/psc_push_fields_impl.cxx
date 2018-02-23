
#include "psc_push_fields_private.h"

#include "psc_fields_single.h"
#include "psc_fields_c.h"

#include "push_fields.hxx"
#include "psc_push_fields_impl.hxx"

// ----------------------------------------------------------------------

template<typename PushFields_t>
static void psc_push_fields_sub_setup(struct psc_push_fields *push)
{
  PscPushFields<PushFields_t> pushf(push);
  new(pushf.sub()) PushFields_t;
}

template<typename PushFields_t>
static void psc_push_fields_sub_destroy(struct psc_push_fields *push)
{
  PscPushFields<PushFields_t> pushf(push);
  pushf.sub()->~PushFields_t();
}

template<typename PushFields_t>
static void psc_push_fields_sub_push_mflds_E(struct psc_push_fields *push,
					     struct psc_mfields *mflds_base,
					     double dt_fac)
{
  PscPushFields<PushFields_t> pushf(push);
  pushf->push_E(push, mflds_base, dt_fac);
}

template<typename PushFields_t>
static void psc_push_fields_sub_push_mflds_H(struct psc_push_fields *push,
					     struct psc_mfields *mflds_base,
					     double dt_fac)
{
  PscPushFields<PushFields_t> pushf(push);
  pushf->push_H(push, mflds_base, dt_fac);
}

// ======================================================================
// psc_push_fields: subclass "single"

struct psc_push_fields_ops_single : psc_push_fields_ops {
  psc_push_fields_ops_single() {
    using PushFields_t = PushFields<mfields_single_t>;
    name                  = "single";
    size                  = sizeof(PushFields_t);
    setup                 = psc_push_fields_sub_setup<PushFields_t>;
    destroy               = psc_push_fields_sub_destroy<PushFields_t>;
    push_mflds_E          = psc_push_fields_sub_push_mflds_E<PushFields_t>;
    push_mflds_H          = psc_push_fields_sub_push_mflds_H<PushFields_t>;
  }
} psc_push_fields_single_ops;

// ======================================================================
// psc_push_fields: subclass "c"

struct psc_push_fields_ops_c : psc_push_fields_ops {
  psc_push_fields_ops_c() {
    using PushFields_t = PushFields<mfields_single_t>;
    name                  = "c";
    size                  = sizeof(PushFields_t);
    setup                 = psc_push_fields_sub_setup<PushFields_t>;
    destroy               = psc_push_fields_sub_destroy<PushFields_t>;
    push_mflds_E          = psc_push_fields_sub_push_mflds_E<PushFields_t>;
    push_mflds_H          = psc_push_fields_sub_push_mflds_H<PushFields_t>;
  }
} psc_push_fields_c_ops;
