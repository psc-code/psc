
#include "vpic/PscRng.h" // FIXME path

#define VPIC_CONFIG_H // FIXME, bad hack

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

#include "rngpool_iface.h"

// ----------------------------------------------------------------------
// RngPool

RngPool* RngPool_create()
{
  return new RngPool;
}

void RngPool_delete(RngPool* rngpool)
{
  delete rngpool;
}

void RngPool_seed(RngPool* rngpool, int base)
{
  rngpool->seed(base, 0);
}

Rng* RngPool_get(RngPool* rngpool, int n)
{
  return (*rngpool)[n];
}

// ----------------------------------------------------------------------
// Rng

double Rng_uniform(Rng *rng, double lo, double hi)
{
  return rng->uniform(lo, hi);
}

double Rng_normal(Rng *rng, double mu, double sigma)
{
  return rng->normal(mu, sigma);
}




