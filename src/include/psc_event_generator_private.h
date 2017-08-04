
#ifndef PSC_EVENT_GENERATOR_PRIVATE_H
#define PSC_EVENT_GENERATOR_PRIVATE_H

#include <psc_event_generator.h>

struct psc_event_generator {
  struct mrc_obj obj;
};

struct psc_event_generator_ops {
  MRC_SUBCLASS_OPS(struct psc_event_generator);
  void (*run)(struct psc_event_generator *event_generator,
	      struct psc_mparticles *mparticles, mfields_base_t *mflds,
	      mphotons_t *mphotons);
};

// ======================================================================

extern struct psc_event_generator_ops psc_event_generator_none_ops;
extern struct psc_event_generator_ops psc_event_generator_demo_ops;

#define psc_event_generator_ops(event_generator) ((struct psc_event_generator_ops *)((event_generator)->obj.ops))

#endif
