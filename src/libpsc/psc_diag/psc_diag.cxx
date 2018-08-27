
#include "psc_diag_private.h"
#include "psc_diag_item_private.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

// ----------------------------------------------------------------------
// psc_diag_setup

static void
_psc_diag_setup(struct psc_diag *diag)
{
  // parse "items" parameter
  char *s_orig = strdup(diag->items), *p, *s = s_orig;
  while ((p = strsep(&s, ", "))) {
    struct psc_diag_item *item =
      psc_diag_item_create(psc_diag_comm(diag));
    psc_diag_item_set_type(item, p);
    psc_diag_item_set_name(item, p);
    psc_diag_add_child(diag, (struct mrc_obj *) item);
  }
  free(s);
  free(s_orig);

  int rank;
  MPI_Comm_rank(psc_diag_comm(diag), &rank);

  if (rank != 0)
    return;
  
  diag->file = fopen("diag.asc", "w");
  fprintf(diag->file, "# time");

  struct psc_diag_item *item;
  mrc_obj_for_each_child(item, diag, struct psc_diag_item) {
    int nr_values = psc_diag_item_nr_values(item);
    for (int i = 0; i < nr_values; i++) {
      fprintf(diag->file, " %s", psc_diag_item_title(item, i));
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
psc_diag_run(struct psc_diag *diag, MparticlesBase& mprts, MfieldsStateBase& mflds)
{
  const auto& grid = mprts.grid();
  
  if (diag->every_step < 0 || 
      grid.timestep() % diag->every_step != 0)
    return;

  int rank;
  MPI_Comm_rank(psc_diag_comm(diag), &rank);

  struct psc_diag_item *item;
  if (rank == 0) {
    fprintf(diag->file, "%g", grid.timestep() * grid.dt);
  }
  mrc_obj_for_each_child(item, diag, struct psc_diag_item) {
    int nr_values = psc_diag_item_nr_values(item);
    double *result = (double *) calloc(nr_values, sizeof(*result));
    psc_diag_item_run(item, mprts, mflds, result);
    if (rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, result, nr_values, MPI_DOUBLE, MPI_SUM, 0, grid.comm());
    } else {
      MPI_Reduce(result, NULL, nr_values, MPI_DOUBLE, MPI_SUM, 0, grid.comm());
    }
    if (rank == 0) {
      for (int i = 0; i < nr_values; i++) {
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
  { "items"            , VAR(items)              ,
    PARAM_STRING("field_energy,particle_energy")                              },
  { "every_step"       , VAR(every_step)         , PARAM_INT(-1)              },
  {},
};
#undef VAR

// ======================================================================
// psc_diag class

struct mrc_class_psc_diag_ : mrc_class_psc_diag {
  mrc_class_psc_diag_() {
    name             = "psc_diag";
    size             = sizeof(struct psc_diag);
    param_descr      = psc_diag_descr;
    setup            = _psc_diag_setup;
    destroy          = _psc_diag_destroy;
  }
} mrc_class_psc_diag;

