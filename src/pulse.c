
#include "psc_pulse.h"

struct psc_pulse *
psc_pulse_create(struct psc_pulse_ops *ops)
{
  struct psc_pulse *pulse = malloc(sizeof(*pulse));
  memset(pulse, 0, sizeof(*pulse));
  pulse->ops = ops;
  
  if (ops->ctx_size) {
    pulse->ctx = malloc(ops->ctx_size);
    memset(pulse->ctx, 0, ops->ctx_size);
  }

  return pulse;
}

void
psc_pulse_destroy(struct psc_pulse *pulse)
{
  if (pulse->ops->destroy) {
    pulse->ops->destroy(pulse);
  }
  free(pulse->ctx);

  free(pulse);
}
