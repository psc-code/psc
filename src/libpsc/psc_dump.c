
#include "psc.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

#include <stdlib.h>
#include <string.h>

// debugging stuff...

static void
ascii_dump_field_yz(struct psc_fields *pf, int m, FILE *file)
{
  for (int iz = pf->ib[2]; iz < pf->ib[2] + pf->im[2]; iz++) {
    for (int iy = pf->ib[1]; iy < pf->ib[1] + pf->im[1]; iy++) {
      int ix = 0; {
	fprintf(file, "%d %d %d %g\n", ix, iy, iz, F3(pf, m, ix,iy,iz));
      }
    }
    fprintf(file, "\n");
  }
  fprintf(file, "\n");
}

static void
ascii_dump_field(struct psc_mfields *flds_base, int m, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (int p = 0; p < flds_base->nr_patches; p++) {
    char *filename = malloc(strlen(fname) + 20);
    sprintf(filename, "%s-p%d-p%d.asc", fname, rank, p);
    mpi_printf(MPI_COMM_WORLD, "ascii_dump_field: '%s'\n", filename);

    struct psc_fields *pf_base = psc_mfields_get_patch(flds_base, p);
    struct psc_fields *pf = psc_fields_get_as(pf_base, "c", m, m+1);
    FILE *file = fopen(filename, "w");
    free(filename);
    if (pf->im[0] + 2*pf->ib[0] == 1) {
      ascii_dump_field_yz(pf, m, file);
    } else {
      for (int iz = pf->ib[2]; iz < pf->ib[2] + pf->im[2]; iz++) {
	for (int iy = pf->ib[1]; iy < pf->ib[1] + pf->im[1]; iy++) {
	  for (int ix = pf->ib[0]; ix < pf->ib[0] + pf->im[0]; ix++) {
	    fprintf(file, "%d %d %d %g\n", ix, iy, iz, F3(pf, m, ix,iy,iz));
	  }
	  fprintf(file, "\n");
	}
	fprintf(file, "\n");
      }
    }
    fclose(file);
    psc_fields_put_as(pf, pf_base, 0, 0);
  }
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
      fprintf(file, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	      i, p->xi, p->yi, p->zi,
	      p->pxi, p->pyi, p->pzi, p->qni, p->mni, p->wni);
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

