
#ifndef RNGPOOL_IFACE_H
#define RNGPOOL_IFACE_H

#include <mrc_common.h>

#ifndef VPIC_CONFIG_H
typedef struct RngPool_ RngPool;
typedef struct Rng_ Rng;
#endif

BEGIN_C_DECLS

RngPool* RngPool_create();
void RngPool_delete(RngPool* rngpool);
void RngPool_seed(RngPool* rngpool, int base);
Rng* RngPool_get(RngPool* rngpool, int n);
double Rng_uniform(Rng *rng, double lo, double hi);
double Rng_normal(Rng *rng, double mu, double sigma);

END_C_DECLS

#endif
