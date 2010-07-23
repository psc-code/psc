#include "psc.h"
#include "output_fields.h"

#include <stdlib.h>
#include <assert.h>

void
init_output_fields(struct psc_extra_fields *f)
{
  f->size = ((psc.ihi[0]-psc.ilo[0]) *
	     (psc.ihi[1]-psc.ilo[1]) *
	     (psc.ihi[2]-psc.ilo[2]));

  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    f->all[m] = malloc(f->size * sizeof(float));
  }

  reset_fields(f);
}

void
free_output_fields(struct psc_extra_fields *f)
{
  for (int m = 0; m < NR_EXTRA_FIELDS ; m++)
    free(f->all[m]);
}

void
calculate_pfields(struct psc_extra_fields *p)
{
   int j = 0;

   for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
      for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
        for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++) {

           p->all[X_EX][j] = .5f * ( FF3(EX,ix,iy,iz)
                                    +FF3(EX,ix-1,iy,iz));
           p->all[X_EY][j] = .5f * ( FF3(EX,ix,iy,iz)
                                    +FF3(EY,ix,iy-1,iz));
           p->all[X_EZ][j] = .5f * ( FF3(EZ,ix,iy,iz)
                                    +FF3(EZ,ix,iy,iz-1));

           p->all[X_BX][j] =  .25f * ( FF3(BX,ix,iy,iz)
                                      +FF3(BX,ix,iy-1,iz)
                                      +FF3(BX,ix,iy,iz-1) 
  			              +FF3(BX,ix,iy-1,iz-1));
           p->all[X_BY][j] =  .25f * ( FF3(BY,ix,iy,iz)
                                      +FF3(BY,ix-1,iy,iz)
                                      +FF3(BY,ix,iy,iz-1) 
				      +FF3(BY,ix-1,iy,iz-1));
           p->all[X_BZ][j] =  .25f * ( FF3(BZ,ix,iy,iz)
                                      +FF3(BZ,ix-1,iy,iz)
                                      +FF3(BZ,ix,iy-1,iz) 
                                      +FF3(BZ,ix-1,iy-1,iz));

           p->all[X_JXI][j] = .5f * ( FF3(JXI,ix,iy,iz) + FF3(JXI,ix-1,iy,iz) );
           p->all[X_JYI][j] = .5f * ( FF3(JYI,ix,iy,iz) + FF3(JYI,ix,iy-1,iz) );
           p->all[X_JZI][j] = .5f * ( FF3(JZI,ix,iy,iz) + FF3(JZI,ix,iy,iz-1) );

           p->all[X_JXEX][j] = p->all[X_JXI][j] * p->all[X_EX][j];
           p->all[X_JYEY][j] = p->all[X_JYI][j] * p->all[X_EY][j];
           p->all[X_JZEZ][j] = p->all[X_JZI][j] * p->all[X_EZ][j];

           p->all[X_POYX][j] = p->all[X_EY][j] * p->all[X_BZ][j] - p->all[X_EZ][j] * p->all[X_BY][j];
           p->all[X_POYY][j] = p->all[X_EZ][j] * p->all[X_BX][j] - p->all[X_EX][j] * p->all[X_BZ][j];
           p->all[X_POYZ][j] = p->all[X_EX][j] * p->all[X_BY][j] - p->all[X_EY][j] * p->all[X_BX][j];

           p->all[X_E2X][j] = p->all[X_EX][j]*p->all[X_EX][j];
           p->all[X_E2Y][j] = p->all[X_EY][j]*p->all[X_EY][j];
           p->all[X_E2Z][j] = p->all[X_EZ][j]*p->all[X_EZ][j];

           p->all[X_B2X][j] = p->all[X_BX][j]*p->all[X_BX][j];
           p->all[X_B2Y][j] = p->all[X_BY][j]*p->all[X_BY][j];
           p->all[X_B2Z][j] = p->all[X_BZ][j]*p->all[X_BZ][j];

           j++;
         }
      }
    }
}

void
accumulate_tfields(struct psc_extra_fields *p, struct psc_extra_fields *t)
{
  t->naccum = t->naccum + 1;

  assert(p->size == t->size);
  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    for (int j = 0; j < p->size; j++)  {
      t->all[m][j] += p->all[m][j];
    }
  }
}


void
reset_fields(struct psc_extra_fields *f)
{
  f->naccum = 0;

  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    for (int j = 0; j < f->size; j++)  {
      f->all[m][j] = 0.0;
    }
  }
}

// convert accumulated values to correct temporal mean
// (divide by naccum)
void
mean_tfields(struct psc_extra_fields *f)
{
  assert(f->naccum > 0);
  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    for (int j = 0; j < f->size; j++) {
      f->all[m][j] /= f->naccum;
    }
  }
}
