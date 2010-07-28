
#include "psc_pulse.h"

struct psc_pulse *
psc_pulse_create(size_t size, struct psc_pulse_ops *ops)
{
  struct psc_pulse *pulse = malloc(sizeof(*pulse));
  memset(pulse, 0, sizeof(*pulse));
  pulse->ops = ops;
  pulse->ctx = malloc(size);
  memset(pulse->ctx, 0, size);

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
