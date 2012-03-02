
#include "psc_diag_private.h"

#include <mrc_params.h>
#include <stdlib.h>

// ----------------------------------------------------------------------
// psc_diag_set_items

void
psc_diag_set_items(struct psc_diag *diag, struct psc_diag_item **items)
{
  diag->items = items;
}

// ----------------------------------------------------------------------
// psc_diag_create

static void
_psc_diag_create(struct psc_diag *diag)
{
  // set default: energy

  static struct psc_diag_item *items[] = {
    &psc_diag_item_em_energy,
    &psc_diag_item_particle_energy,
    NULL,
  };
  
  diag->items = items;
}

// ----------------------------------------------------------------------
// psc_diag_setup

static void
_psc_diag_setup(struct psc_diag *diag)
{
  int rank;
  MPI_Comm_rank(psc_diag_comm(diag), &rank);

  if (rank != 0)
    return;
  
  diag->file = fopen("diag.asc", "w");
  fprintf(diag->file, "# time");
  for (int m = 0; diag->items[m]; m++) {
    struct psc_diag_item *item = diag->items[m];
    for (int i = 0; i < item->n_values; i++) {
      fprintf(diag->file, " %s", item->names[i]);
    }
  }
  fprintf(diag->file, "\n");
}

// ----------------------------------------------------------------------
// psc_diag_destroy

static void
_psc_diag_destroy(struct psc_diag *diag)
{
  int rank;
  MPI_Comm_rank(psc_diag_comm(diag), &rank);

  if (rank != 0)
    return;
  
  fclose(diag->file);
}

// ----------------------------------------------------------------------
// psc_diag_run

void
psc_diag_run(struct psc_diag *diag, struct psc *psc)
{
  if (diag->every_step < 0 || 
      psc->timestep % diag->every_step != 0)
    return;

  int rank;
  MPI_Comm_rank(psc_diag_comm(diag), &rank);

  for (int m = 0; diag->items[m]; m++) {
    struct psc_diag_item *item = diag->items[m];

    double *result = calloc(item->n_values, sizeof(*result));
    item->run(psc, result);
    MPI_Reduce(MPI_IN_PLACE, result, item->n_values, MPI_DOUBLE, MPI_SUM, 0, psc_comm(psc));
    if (rank == 0) {
      fprintf(diag->file, "%g", psc->timestep * psc->dt);
      for (int i = 0; i < item->n_values; i++) {
	fprintf(diag->file, " %g", result[i]);
      }
    }
    free(result);
  }
  if (rank == 0) {
    fprintf(diag->file, "\n");
    fflush(diag->file);
  }
}

// ======================================================================

#define VAR(x) (void *)offsetof(struct psc_diag, x)

static struct param psc_diag_descr[] = {
  { "every_step"       , VAR(every_step)         , PARAM_INT(-1)            },
  {},
};
#undef VAR

// ======================================================================
// psc_diag class

struct mrc_class_psc_diag mrc_class_psc_diag = {
  .name             = "psc_diag",
  .size             = sizeof(struct psc_diag),
  .param_descr      = psc_diag_descr,
  .create           = _psc_diag_create,
  .setup            = _psc_diag_setup,
  .destroy          = _psc_diag_destroy,
};

