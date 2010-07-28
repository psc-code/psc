
#include "psc_pulse.h"

struct psc_pulse *
psc_pulse_create(struct psc_pulse_ops *ops, void *prm)
{
  struct psc_pulse *pulse = malloc(sizeof(*pulse));
  memset(pulse, 0, sizeof(*pulse));
  pulse->ops = ops;
  
  if (ops->ctx_size) {
    pulse->ctx = malloc(ops->ctx_size);
    memset(pulse->ctx, 0, ops->ctx_size);
  }

  struct param *descr = ops->ctx_descr;
  if (descr) {
    char name[strlen(ops->name) + 7];
    sprintf(name, "Pulse %s", ops->name);
    if (prm) { // custom defaults were passed
      memcpy(pulse->ctx, prm, ops->ctx_size);
      params_parse_cmdline_nodefault(pulse->ctx, descr, name, MPI_COMM_WORLD);
    } else {
      params_parse_cmdline(pulse->ctx, descr, name, MPI_COMM_WORLD);
    }
    params_print(pulse->ctx, descr, name, MPI_COMM_WORLD);
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
