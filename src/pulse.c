
#include "psc_pulse.h"

struct psc_pulse *
psc_pulse_create(size_t size, struct psc_pulse_ops *ops)
{
  struct psc_pulse *pulse = malloc(size);
  memset(pulse, 0, size);
  pulse->ops = ops;

  return pulse;
}

void
psc_pulse_destroy(struct psc_pulse *pulse)
{
  if (pulse->ops->destroy) {
    pulse->ops->destroy(pulse);
  }

  free(pulse);
}
