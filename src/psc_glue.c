
#include "psc.h"

#define PIC_push_part_yz_F77 F77_FUNC(pic_push_part_yz,PIC_PUSH_PART_YZ)

void PIC_push_part_yz_F77(f_int *niloc, struct f_particle *p_niloc,
			  f_real *cori, f_real *alpha, f_real *eta, 
			  f_real *p2A, f_real *p2B,
			  f_real *dt, f_real *dx, f_real *dy, f_real *dz,
			  f_int *i1mn, f_int *i2mn, f_int *i3mn,
			  f_int *i1mx, f_int *i2mx, f_int *i3mx,
			  f_int *rd1, f_int *rd2, f_int *rd3,
			  f_real *ne, f_real *ni, f_real *nn,
			  f_real *jxi, f_real *jyi, f_real *jzi,
			  f_real *ex, f_real *ey, f_real *ez,
			  f_real *bx, f_real *by, f_real *bz);

void
PIC_push_part_yz()
{
  int i0mx = psc.ihi[0] - 1, i1mx = psc.ihi[1] - 1, i2mx = psc.ihi[2] - 1;

  PIC_push_part_yz_F77(&psc.n_part, &psc.f_part[-1], &psc.prm.cori, &psc.prm.alpha,
		       &psc.prm.eta, &psc.p2A, &psc.p2B,
		       &psc.dt, &psc.dx[0], &psc.dx[1], &psc.dx[2],
		       &psc.ilo[0], &psc.ilo[1], &psc.ilo[2],
		       &i0mx, &i1mx, &i2mx,
		       &psc.ibn[0], &psc.ibn[1], &psc.ibn[2],
		       psc.f_fields[NE], psc.f_fields[NI], psc.f_fields[NN],
		       psc.f_fields[JXI], psc.f_fields[JYI], psc.f_fields[JZI],
		       psc.f_fields[EX], psc.f_fields[EY], psc.f_fields[EZ],
		       psc.f_fields[BX], psc.f_fields[BY], psc.f_fields[BZ]);
}

