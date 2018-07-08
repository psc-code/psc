
#pragma once

#include "psc.h"
#include "particles.hxx"

#include "psc_stats.h"
#include <mrc_profile.h>

// ======================================================================
// SortBase

struct SortBase
{
  virtual void run(MparticlesBase& mprts_base) = 0;
};

// ======================================================================
// SortConvert

template<typename Sort_t>
struct SortConvert : Sort_t, SortBase
{
  using Base = Sort_t;
  using Mparticles = typename Sort_t::Mparticles;
  using Base::Base;

  void run(MparticlesBase& mprts_base) override
  {
    auto& mprts = mprts_base.get_as<Mparticles>();
    (*this)(mprts);
    mprts_base.put_as(mprts);
  }
};

