
#include "vpic/PscRng.h" // FIXME path

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

#include "rngpool_iface.h"


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



