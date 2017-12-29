
#ifndef RNGPOOL_IFACE_H
#define RNGPOOL_IFACE_H

#include <mrc_common.h>

#ifndef __cplusplus
typedef struct RngPool_ RngPool;
typedef struct Rng_ Rng;
#endif

BEGIN_C_DECLS

RngPool* RngPool_create();
void RngPool_delete(RngPool* rngpool);
void RngPool_seed(RngPool* rngpool, int base);
Rng* RngPool_get(RngPool* rngpool, int n);

END_C_DECLS

#endif
