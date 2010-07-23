#include "psc.h"
#include "output_fields.h"

#include <stdlib.h>
#include <assert.h>

void
init_output_fields(struct psc_extra_fields *f)
{
  f->size =  (psc.ihi[0]-psc.ilo[0])
                      * (psc.ihi[1]-psc.ilo[1])
                      * (psc.ihi[2]-psc.ilo[2]);
  f->nfields = 21; 

  f->ex = f->all[0] = malloc(f->size * sizeof(float) );
  f->ey = f->all[1] = malloc(f->size * sizeof(float) );
  f->ez = f->all[2] = malloc(f->size * sizeof(float) );
  f->bx = f->all[3] = malloc(f->size * sizeof(float) );
  f->by = f->all[4] = malloc(f->size * sizeof(float) );
  f->bz = f->all[5] = malloc(f->size * sizeof(float) );
  f->jx = f->all[6] = malloc(f->size * sizeof(float) );
  f->jy = f->all[7] = malloc(f->size * sizeof(float) );
  f->jz = f->all[8] = malloc(f->size * sizeof(float) );
  f->jxex = f->all[9]  = malloc(f->size * sizeof(float) );
  f->jyey = f->all[10] = malloc(f->size * sizeof(float) );
  f->jzez = f->all[11] = malloc(f->size * sizeof(float) );
  f->poyx = f->all[12] = malloc(f->size * sizeof(float) );
  f->poyy = f->all[13] = malloc(f->size * sizeof(float) );
  f->poyz = f->all[14] = malloc(f->size * sizeof(float) );
  f->e2x = f->all[15] = malloc(f->size * sizeof(float) );
  f->e2y = f->all[16] = malloc(f->size * sizeof(float) );
  f->e2z = f->all[17] = malloc(f->size * sizeof(float) );
  f->b2x = f->all[18] = malloc(f->size * sizeof(float) );
  f->b2y = f->all[19] = malloc(f->size * sizeof(float) );
  f->b2z = f->all[20] = malloc(f->size * sizeof(float) );

  reset_fields(f);
}

void
free_output_fields(struct psc_extra_fields *f)
{
  for (int m = 0; m<f->nfields ; m++)
    free(f->all[m]);
}

void
calculate_pfields(struct psc_extra_fields *p)
{
   int j = 0;

   for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
      for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
        for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++) {

           p->ex[j] = (float) 0.5  * ( FF3(EX,ix,iy,iz)
                                     +FF3(EX,ix-1,iy,iz) );
           p->ey[j] = (float) 0.5  * ( FF3(EX,ix,iy,iz)
                                     +FF3(EY,ix,iy-1,iz) );
           p->ez[j] = (float) 0.5  * ( FF3(EZ,ix,iy,iz)
                                     +FF3(EZ,ix,iy,iz-1) );

           p->bx[j] = (float) 0.25 * ( FF3(BX,ix,iy,iz)
                                     +FF3(BX,ix,iy-1,iz)
                                     +FF3(BX,ix,iy,iz-1) 
                                     +FF3(BX,ix,iy-1,iz-1) );
           p->by[j] = (float) 0.25 * ( FF3(BY,ix,iy,iz)
                                     +FF3(BY,ix-1,iy,iz)
                                     +FF3(BY,ix,iy,iz-1) 
                                     +FF3(BY,ix-1,iy,iz-1) );
           p->bz[j] = (float) 0.25 * ( FF3(BZ,ix,iy,iz)
                                     +FF3(BZ,ix-1,iy,iz)
                                     +FF3(BZ,ix,iy-1,iz) 
                                     +FF3(BZ,ix-1,iy-1,iz) );

           p->jx[j] = (float) 0.5  * ( FF3(JXI,ix,iy,iz)
                                     +FF3(JXI,ix-1,iy,iz) );
           p->jy[j] = (float) 0.5  * ( FF3(JYI,ix,iy,iz)
                                     +FF3(JYI,ix,iy-1,iz) );
           p->jz[j] = (float) 0.5  * ( FF3(JZI,ix,iy,iz)
                                     +FF3(JZI,ix,iy,iz-1) );

           p->jxex[j] = p->jx[j] * p->ex[j];
           p->jyey[j] = p->jy[j] * p->ey[j];
           p->jzez[j] = p->jz[j] * p->ez[j];

           p->poyx[j] = p->ey[j] * p->bz[j] - p->ez[j]*p->by[j];
           p->poyy[j] = p->ez[j] * p->bx[j] - p->ex[j]*p->bz[j];
           p->poyz[j] = p->ex[j] * p->by[j] - p->ey[j]*p->bx[j];

           p->e2x[j] = p->ex[j]*p->ex[j];
           p->e2y[j] = p->ey[j]*p->ey[j];
           p->e2z[j] = p->ez[j]*p->ez[j];

           p->b2x[j] = p->bx[j]*p->bx[j];
           p->b2y[j] = p->by[j]*p->by[j];
           p->b2z[j] = p->bz[j]*p->bz[j];

           j=j+1;
         }
      }
    }
}

void
accumulate_tfields(struct psc_extra_fields *p, struct psc_extra_fields *t)
{
  t->naccum = t->naccum + 1;

  assert(p->size == t->size);
  for (int m=0; m < p->nfields; m++) {
    for (int j=0; j < p->size; j++)  {
      t->all[m][j] += p->all[m][j];
    }
  }
}


void
reset_fields(struct psc_extra_fields *f)
{
  f->naccum = 0;

  for (int m=0; m < f->nfields; m++) {
    for (int j=0; j < f->size; j++)  {
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
  for (int m=0; m < f->nfields; m++) {
    for (int j=0; j  <f->size; j++) {
      f->all[m][j] = f->all[m][j] / f->naccum;
    }
  }
}
