
#pragma once

#include "psc_output_fields_private.h"

// ======================================================================
// PscOutputFields

template<typename S>
struct PscOutputFields
{
  using sub_t = S;

  explicit PscOutputFields(psc_output_fields* outf)
    : outf_{outf}
  {}

  sub_t* sub() { return mrc_to_subobj(outf_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_output_fields* outf_;
};
