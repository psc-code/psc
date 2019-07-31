
#include "psc_diag_item_private.h"
#include "psc_diag_private.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

// ----------------------------------------------------------------------
// psc_diag_setup

static void _psc_diag_setup(struct psc_diag* diag)
{
  struct psc_diag_item* item = psc_diag_item_create(psc_diag_comm(diag));
  psc_diag_item_set_type(item, "field_energy");
  psc_diag_item_set_name(item, "field_energy");
  diag->items_.push_back(item);

  item = psc_diag_item_create(psc_diag_comm(diag));
  psc_diag_item_set_type(item, "particle_energy");
  psc_diag_item_set_name(item, "particle_energy");
  diag->items_.push_back(item);

  MPI_Comm_rank(psc_diag_comm(diag), &diag->rank_);

  if (diag->rank_ == 0) {
    diag->file_ = fopen("diag.asc", "w");
    fprintf(diag->file_, "# time");
    
    for (auto item: diag->items_) {
      int nr_values = psc_diag_item_nr_values(item);
      for (int i = 0; i < nr_values; i++) {
	fprintf(diag->file_, " %s", psc_diag_item_title(item, i));
      }
    }
    fprintf(diag->file_, "\n");
  }
}

// ----------------------------------------------------------------------
// psc_diag_destroy

static void _psc_diag_destroy(struct psc_diag* diag)
{
  if (diag->rank_ == 0) {
    fclose(diag->file_);
  }
}

// ----------------------------------------------------------------------
// psc_diag_run

void psc_diag_run(struct psc_diag* diag, MparticlesBase& mprts,
                  MfieldsStateBase& mflds)
{
  const auto& grid = mprts.grid();

  if (diag->every_step < 0 || grid.timestep() % diag->every_step != 0)
    return;

  if (diag->rank_ == 0) {
    fprintf(diag->file_, "%g", grid.timestep() * grid.dt);
  }
  for (auto item: diag->items_) {
    int nr_values = psc_diag_item_nr_values(item);
    std::vector<double> result(nr_values);
    psc_diag_item_run(item, mprts, mflds, result.data());
    if (diag->rank_ == 0) {
      MPI_Reduce(MPI_IN_PLACE, result.data(), result.size(), MPI_DOUBLE, MPI_SUM, 0,
                 grid.comm());
    } else {
      MPI_Reduce(result.data(), NULL, result.size(), MPI_DOUBLE, MPI_SUM, 0, grid.comm());
    }
    if (diag->rank_ == 0) {
      for (auto val : result) {
        fprintf(diag->file_, " %g", val);
      }
    }
  }
  if (diag->rank_ == 0) {
    fprintf(diag->file_, "\n");
    fflush(diag->file_);
  }
}

// ======================================================================

#define VAR(x) (void*)offsetof(struct psc_diag, x)

static struct param psc_diag_descr[] = {
  {"every_step", VAR(every_step), PARAM_INT(-1)},
  {},
};
#undef VAR

// ======================================================================
// psc_diag class

struct mrc_class_psc_diag_ : mrc_class_psc_diag
{
  mrc_class_psc_diag_()
  {
    name = "psc_diag";
    size = sizeof(struct psc_diag);
    param_descr = psc_diag_descr;
    setup = _psc_diag_setup;
    destroy = _psc_diag_destroy;
  }
} mrc_class_psc_diag;
