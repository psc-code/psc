
#include "psc_event_generator_private.h"

// ======================================================================
// forward to subclass

void
psc_event_generator_run(struct psc_event_generator *gen,
			struct psc_mparticles *mparticles, mfields_base_t *mflds)
{
  struct psc_event_generator_ops *ops = psc_event_generator_ops(gen);
  assert(ops->run);
  ops->run(gen, mparticles, mflds);
}

// ======================================================================
// psc_event_generator_init

static void
psc_event_generator_init()
{
  mrc_class_register_subclass(&mrc_class_psc_event_generator, &psc_event_generator_none_ops);
}

// ======================================================================
// psc_event_generator class

struct mrc_class_psc_event_generator mrc_class_psc_event_generator = {
  .name             = "psc_event_generator",
  .size             = sizeof(struct psc_event_generator),
  .init             = psc_event_generator_init,
};

