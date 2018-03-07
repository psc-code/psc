
#pragma once

#include "psc_bnd_private.h"
#include "fields3d.hxx"

#include <mrc_profile.h>
#include <mrc_ddc.h>

// ======================================================================
// BndBase

struct BndBase
{
  virtual void fill_ghosts(PscMfieldsBase mflds, int mb, int me) = 0;
  virtual void add_ghosts(PscMfieldsBase mflds, int mb, int me) = 0;
  virtual void reset() = 0;
};

// ======================================================================
// PscBnd

template<typename S>
struct PscBnd
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, BndBase*>::value,
  		"sub classes used in PscBndParticles must derive from BndBase");
  
  explicit PscBnd(psc_bnd *bnd)
    : bnd_(bnd)
  {}

  void fill_ghosts(PscMfieldsBase mflds, int mb, int me)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("fill_ghosts", 1., 0, 0);
    }
    
    psc_stats_start(st_time_comm);
    prof_start(pr);
    
    sub()->fill_ghosts(mflds.mflds(), mb, me);

    prof_stop(pr);
    psc_stats_stop(st_time_comm);
  }

  void add_ghosts(PscMfieldsBase mflds, int mb, int me)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("add_ghosts", 1., 0, 0);
    }

    psc_stats_start(st_time_comm);
    prof_start(pr);
    
    sub()->add_ghosts(mflds.mflds(), mb, me);
    
    prof_stop(pr);
    psc_stats_stop(st_time_comm);
  }

  void reset()
  {
    sub()->reset();
  }
  
  sub_t* sub() { return mrc_to_subobj(bnd_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_bnd* bnd_;
};

using PscBndBase = PscBnd<BndBase>;

// ======================================================================
// PscBndWrapper

template<typename Bnd>
class PscBndWrapper
{
public:
  const static size_t size = sizeof(Bnd);
  
  static void setup(psc_bnd* _bnd)
  {
    PscBnd<Bnd> bnd(_bnd);
    new(bnd.sub()) Bnd(_bnd->psc->grid(), _bnd->psc->mrc_domain, _bnd->psc->ibn);
  }

  static void destroy(psc_bnd* _bnd)
  {
    PscBnd<Bnd> bnd(_bnd);
    bnd->~Bnd();
  }
};

