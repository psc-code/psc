
#pragma once

#include <psc_method.h>

struct Psc
{
  // ----------------------------------------------------------------------
  // initialize

  void initialize(psc* psc_)
  {
    psc_method_initialize(psc_->method, psc_);
  }
};
