
#include "psc_bnd_private.h"

#include <mrc_io.h>

// ----------------------------------------------------------------------
// psc_set_psc

void
psc_bnd_set_psc(struct psc_bnd *bnd, struct psc *psc)
{
  bnd->psc = psc;
}

// ----------------------------------------------------------------------
// psc_bnd_write

static void
_psc_bnd_write(struct psc_bnd *bnd, struct mrc_io *io)
{
  const char *path = psc_bnd_name(bnd);
  mrc_io_write_obj_ref(io, path, "psc", (struct mrc_obj *) bnd->psc);
}

// ----------------------------------------------------------------------
// psc_bnd_read

static void
_psc_bnd_read(struct psc_bnd *bnd, struct mrc_io *io)
{
  const char *path = psc_bnd_name(bnd);
  bnd->psc = (struct psc *)
    mrc_io_read_obj_ref(io, path, "psc", &mrc_class_psc);
}

// ======================================================================
// forward to subclass

void
psc_bnd_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me)
{
  psc_stats_start(st_time_comm);
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->add_ghosts);
  ops->add_ghosts(bnd, flds, mb, me);
  psc_stats_stop(st_time_comm);
}

void
psc_bnd_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me)
{
  psc_stats_start(st_time_comm);
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->fill_ghosts);
  ops->fill_ghosts(bnd, flds, mb, me);
  psc_stats_stop(st_time_comm);
}

void
psc_bnd_exchange_particles(struct psc_bnd *bnd, mparticles_base_t *particles)
{
  psc_stats_start(st_time_comm);
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->exchange_particles);
  ops->exchange_particles(bnd, particles);
  psc_stats_stop(st_time_comm);
}

void
psc_bnd_exchange_photons(struct psc_bnd *bnd, mphotons_t *mphotons)
{
  psc_stats_start(st_time_comm);
  int n_total = 0;
  psc_foreach_patch(bnd->psc, p) {
    n_total += mphotons->p[p].nr;
  }
  MPI_Allreduce(MPI_IN_PLACE, &n_total, 1, MPI_INT, MPI_SUM, 
		psc_bnd_comm(bnd));
  if (n_total == 0)
    return;

  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->exchange_photons);
  ops->exchange_photons(bnd, mphotons);
  psc_stats_stop(st_time_comm);
}

// ======================================================================
// psc_bnd_init

static void
psc_bnd_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_fortran_ops);
}

// ======================================================================
// psc_bnd class

struct mrc_class_psc_bnd mrc_class_psc_bnd = {
  .name             = "psc_bnd",
  .size             = sizeof(struct psc_bnd),
  .init             = psc_bnd_init,
  .write            = _psc_bnd_write,
  .read             = _psc_bnd_read,
};

