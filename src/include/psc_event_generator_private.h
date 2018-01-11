
#ifndef PSC_EVENT_GENERATOR_PRIVATE_H
#define PSC_EVENT_GENERATOR_PRIVATE_H

#include <psc_event_generator.h>

struct psc_event_generator {
  struct mrc_obj obj;
};

struct psc_event_generator_ops {
  MRC_SUBCLASS_OPS(struct psc_event_generator);
  void (*run)(struct psc_event_generator *event_generator,
	      struct psc_mparticles *mprts, struct psc_mfields *mflds);
};

// ======================================================================

#define psc_event_generator_ops(event_generator) ((struct psc_event_generator_ops *)((event_generator)->obj.ops))

#endif
