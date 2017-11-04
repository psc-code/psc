
#include "psc.h"
#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include <stdlib.h>
#include <string.h>

// debugging stuff...

static void
ascii_dump_field_yz(fields_t flds, int m, FILE *file)
{
  for (int iz = flds.ib[2]; iz < flds.ib[2] + flds.im[2]; iz++) {
    for (int iy = flds.ib[1]; iy < flds.ib[1] + flds.im[1]; iy++) {
      int ix = 0; {
	fprintf(file, "%d %d %d %g\n", ix, iy, iz, _F3(flds, m, ix,iy,iz));
      }
    }
    fprintf(file, "\n");
  }
  fprintf(file, "\n");
}

static void
ascii_dump_field(struct psc_mfields *mflds_base, int m, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "c", m, m+1);
  for (int p = 0; p < mflds_base->nr_patches; p++) {
    char *filename = malloc(strlen(fname) + 20);
    sprintf(filename, "%s-p%d-p%d.asc", fname, rank, p);
    mpi_printf(MPI_COMM_WORLD, "ascii_dump_field: '%s'\n", filename);

    fields_t flds = fields_t_mflds(mflds, p);
    FILE *file = fopen(filename, "w");
    free(filename);
    if (flds.im[0] + 2*flds.ib[0] == 1) {
      ascii_dump_field_yz(flds, m, file);
    } else {
      for (int iz = flds.ib[2]; iz < flds.ib[2] + flds.im[2]; iz++) {
	for (int iy = flds.ib[1]; iy < flds.ib[1] + flds.im[1]; iy++) {
	  for (int ix = flds.ib[0]; ix < flds.ib[0] + flds.im[0]; ix++) {
	    fprintf(file, "%d %d %d %g\n", ix, iy, iz, _F3(flds, m, ix,iy,iz));
	  }
	  fprintf(file, "\n");
	}
	fprintf(file, "\n");
      }
    }
    fclose(file);
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

static void
ascii_dump_particles(struct psc_mparticles *mprts_base, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "c", 0);
  
  psc_foreach_patch(ppsc, p) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    char *filename = malloc(strlen(fname) + 20);
    sprintf(filename, "%s-p%d-p%d.asc", fname, rank, p);
    mpi_printf(MPI_COMM_WORLD, "ascii_dump_particles: '%s'\n", filename);
    
    FILE *file = fopen(filename, "w");
    fprintf(file, "i\txi\tyi\tzi\tpxi\tpyi\tpzi\tqni\tmni\twni\n");
    for (int i = 0; i < particle_range_size(prts); i++) {
      particle_t *p = particle_iter_at(prts.begin, i);
      fprintf(file, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n",
	      i, p->xi, p->yi, p->zi,
	      p->pxi, p->pyi, p->pzi, p->qni_wni, p->kind);
    }
    fclose(file);
    free(filename);
  }

  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
}

void
psc_dump_field(struct psc_mfields *flds, int m, const char *fname)
{
  ascii_dump_field(flds, m, fname);
}

void
psc_dump_particles(struct psc_mparticles *particles, const char *fname)
{
  ascii_dump_particles(particles, fname);
}

