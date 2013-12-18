#ifndef MRC_BLOCK_FACTORY_H
#define MRC_BLOCK_FACTORY_H

#include <mrc_domain.h>

MRC_CLASS_DECLARE(mrc_block_factory, struct mrc_block_factory);


void mrc_block_factory_run(struct mrc_block_factory *fac, struct mrc_domain *domain);


#endif
