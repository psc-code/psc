
#include "psc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct psc psc;

void
psc_alloc(int ilo[3], int ihi[3], int ibn[3], int n_part)
{
  for (int d = 0; d < 3; d++) {
    psc.ilo[d] = ilo[d];
    psc.ilg[d] = ilo[d] - ibn[d];
    psc.ihi[d] = ihi[d];
    psc.ihg[d] = ihi[d] + ibn[d];
    psc.img[d] = ihi[d] - ilo[d] + 2 * ibn[d];
  }
  psc.fld_size = psc.img[0] * psc.img[1] * psc.img[2];
  for (int m = 0; m < NR_FIELDS; m++) {
    psc.f_fields[m] = calloc(psc.fld_size, sizeof(f_real));
  }

  psc.n_part = n_part;
  psc.f_part = calloc(n_part, sizeof(*psc.f_part));
}

void
psc_free()
{
  for (int m = 0; m < NR_FIELDS; m++) {
    free(psc.f_fields[m]);
  }

  free(psc.f_part);
}

void
psc_setup_parameters()
{
  psc.prm.cori = 2.;
  psc.prm.eta = 3.;
  psc.prm.alpha = 5.;
  psc.dt = 1.;
  psc.dx[0] = 1.;
  psc.dx[1] = 1.;
  psc.dx[2] = 1.;
}

void
psc_setup_fields_zero()
{
  for (int m = 0; m < NR_FIELDS; m++) {
    memset(psc.f_fields[m], 0, psc.fld_size * sizeof(f_real));
  }
}

void
psc_setup_particles_1()
{
  for (int n = 0; n < psc.n_part; n++) {
    psc.f_part[n].xi = .5;
    psc.f_part[n].yi = .5;
    psc.f_part[n].zi = .5 + n;
    psc.f_part[n].pxi = 0.;
    psc.f_part[n].pyi = .02;
    psc.f_part[n].pzi = .01;
    psc.f_part[n].qni = 0.;
    psc.f_part[n].mni = 1.;
    psc.f_part[n].lni = 0.;
    psc.f_part[n].wni = 1.;
  }
}

void
psc_dump_particles(const char *fname)
{
  printf("psc_dump_particles %s\n", fname);

  FILE *file = fopen(fname, "w");
  fprintf(file, "i\txi\tyi\tzi\tpxi\tpyi\tpzi\tqni\tmni\twni\n");
  for (int i = 0; i < psc.n_part; i++) {
    struct f_particle *p = &psc.f_part[i];
    fprintf(file, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	    i, p->xi, p->yi, p->zi,
	    p->pxi, p->pyi, p->pzi, p->qni, p->mni, p->wni);
  }
  fclose(file);
}
