
#include "psc.h"

#include <stdlib.h>
#include <string.h>

// debugging stuff...

static void
ascii_dump_field(mfields_base_t *flds, int m, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  psc_foreach_patch(ppsc, p) {
    char *filename = malloc(strlen(fname) + 10);
    sprintf(filename, "%s-p%d-p%d.asc", fname, rank, p);
    mpi_printf(MPI_COMM_WORLD, "ascii_dump_field: '%s'\n", filename);

    fields_base_t *pf = &flds->f[p];
    FILE *file = fopen(filename, "w");
    free(filename);
    psc_foreach_patch(ppsc, patch) {
      for (int iz = -ppsc->ibn[2]; iz < ppsc->patch[patch].ldims[2] + ppsc->ibn[2]; iz++) {
	for (int iy = -ppsc->ibn[1]; iy < ppsc->patch[patch].ldims[1] + ppsc->ibn[1]; iy++) {
	  for (int ix = -ppsc->ibn[0]; ix < ppsc->patch[patch].ldims[0] +  ppsc->ibn[0]; ix++) {
	    fprintf(file, "%d %d %d %g\n", ix, iy, iz, F3_BASE(pf, m, ix,iy,iz));
	  }
	  fprintf(file, "\n");
	}
	fprintf(file, "\n");
      }
    }
    fclose(file);
  }
}

static void
ascii_dump_particles(mparticles_base_t *particles, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  psc_foreach_patch(ppsc, p) {
    particles_base_t *pp = &particles->p[p];
    char *filename = malloc(strlen(fname) + 10);
    sprintf(filename, "%s-p%d-p%d.asc", fname, rank, p);
    mpi_printf(MPI_COMM_WORLD, "ascii_dump_particles: '%s'\n", filename);
    
    FILE *file = fopen(filename, "w");
    fprintf(file, "i\txi\tyi\tzi\tpxi\tpyi\tpzi\tqni\tmni\twni\n");
    for (int i = 0; i < pp->n_part; i++) {
      particle_base_t *p = particles_base_get_one(pp, i);
      fprintf(file, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	      i, p->xi, p->yi, p->zi,
	      p->pxi, p->pyi, p->pzi, p->qni, p->mni, p->wni);
    }
    fclose(file);
    free(filename);
  }
}

void
psc_dump_field(mfields_base_t *flds, int m, const char *fname)
{
  ascii_dump_field(flds, m, fname);
}

void
psc_dump_particles(mparticles_base_t *particles, const char *fname)
{
  ascii_dump_particles(particles, fname);
}

