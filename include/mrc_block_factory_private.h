#ifndef MRC_BLOCK_FACTORY_PRIVATE_H
#define MRC_BLOCK_FACTORY_PRIVATE_H

#include <mrc_block_factory.h>
#include <mrc_domain_private.h>

struct mrc_block_factory {
  struct mrc_obj obj;
};
  
struct mrc_block_factory_ops {
  MRC_SUBCLASS_OPS(struct mrc_block_factory);
  void (*run)(struct mrc_block_factory *fac, struct mrc_domain *domain);
			  
};

#define mrc_block_factory_ops(fac) ((struct mrc_block_factory_ops *)(fac)->obj.ops)

extern struct mrc_block_factory_ops mrc_block_factory_simple2d;
extern struct mrc_block_factory_ops mrc_block_factory_simple3d;
extern struct mrc_block_factory_ops mrc_block_factory_cylindrical;
extern struct mrc_block_factory_ops mrc_block_factory_half_cylinder;

#endif
