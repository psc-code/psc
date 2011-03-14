
#include "psc.h"
#include "psc_glue.h"

// FIXME, this all needs to be redone

// ----------------------------------------------------------------------
// psc_write_checkpoint

void
psc_read_checkpoint(void)
{
  mpi_printf(MPI_COMM_WORLD, "INFO: Reading checkpoint.\n");
  
  int n_part;
  SERV_read_1(&psc.timestep, &n_part);
  particles_base_realloc(&psc.particles.p[0], n_part);
  psc.particles.p[0].n_part = n_part;
  
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, &psc.particles);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, 0, 0, &psc.flds);

  SERV_read_2(&particles.p[0], &flds.f[0]);

  particles_fortran_put(&particles, &psc.particles);
  fields_fortran_put(&flds, JXI, HZ + 1, &psc.flds);
}

// ----------------------------------------------------------------------
// psc_write_checkpoint

void
psc_write_checkpoint(void)
{
  mpi_printf(MPI_COMM_WORLD, "INFO: Writing checkpoint.\n");
  
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, &psc.particles);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, JXI, HZ + 1, &psc.flds);

  SERV_write(&particles.p[0], &flds.f[0]);

  particles_fortran_put(&particles, &psc.particles);
  fields_fortran_put(&flds, 0, 0, &psc.flds);
}

