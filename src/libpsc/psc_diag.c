
#include "psc.h"

// ----------------------------------------------------------------------
// psc_diag_run

void
psc_diag_run(struct psc *psc, struct psc_diag_item **diag_items)
{
  static FILE *file;
  int rank;

  MPI_Comm_rank(psc_comm(psc), &rank);

  if (!file && rank == 0) {
    file = fopen("diag.asc", "w");
    fprintf(file, "# time");
    for (int m = 0; diag_items[m]; m++) {
      struct psc_diag_item *item = diag_items[m];
      for (int i = 0; i < item->n_values; i++) {
	fprintf(file, " %s", item->names[i]);
      }
    }
    fprintf(file, "\n");
  }

  for (int m = 0; diag_items[m]; m++) {
    struct psc_diag_item *item = diag_items[m];

    double *result = calloc(item->n_values, sizeof(*result));
    item->run(psc, result);
    MPI_Reduce(MPI_IN_PLACE, result, item->n_values, MPI_DOUBLE, MPI_SUM, 0, psc_comm(psc));
    if (rank == 0) {
      fprintf(file, "%g", psc->timestep * psc->dt);
      for (int i = 0; i < item->n_values; i++) {
	fprintf(file, " %g", result[i]);
      }
    }
    free(result);
  }
  fprintf(file, "\n");
  fflush(file);
}

