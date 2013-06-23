
#include <mrc_nr.h>

#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* Taken from NR */
void
mrc_nr_tridag(float *a, float *b, float *c, float *r, float *u, int n)
{
  int j;
  float *gam, bet;
  
  gam = calloc(sizeof(*gam), n);
  assert(fabs(b[0]) >= 1e-7);
  bet = b[0];
  u[0] = r[0] / bet;
  for (j = 1; j < n; j++) {
    gam[j] = c[j-1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (j < n-1 || fabs(bet) >= 1e-7) {
      assert(fabs(bet) >= 1e-7);
      u[j] = (r[j] - a[j] * u[j-1]) / bet;
    } else {
      assert(fabs(r[j] - a[j] * u[j-1]) < 1e-7);
      u[j] = 0.;
    }
  }
    
  for (j = n - 2; j >= 0; j--) {
    u[j] = u[j] - gam[j+1] * u[j+1];
  }
  free(gam);
}

