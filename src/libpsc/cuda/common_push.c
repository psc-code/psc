
// ----------------------------------------------------------------------
// push_xi_z
//
// advance position using velocity

#if DIM == DIM_Z

__device__ static void
push_xi(struct d_particle *p, const real vxi[3], real dt)
{
  for (int d = 2; d < 3; d++) {
    p->xi[d] += dt * vxi[d];
  }
}

#elif DIM == DIM_YZ

__device__ static void
push_xi(struct d_particle *p, const real vxi[3], real dt)
{
  int d;
  for (d = 1; d < 3; d++) {
    p->xi[d] += dt * vxi[d];
  }
}

#else
#error unknown DIM
#endif

// ----------------------------------------------------------------------
// calc_vxi
//
// calculate velocity from moments

__device__ static void
calc_vxi(real vxi[3], struct d_particle p)
{
  real root = rsqrtr(real(1.) + sqr(p.pxi[0]) + sqr(p.pxi[1]) + sqr(p.pxi[2]));

  int d;
  for (d = 0; d < 3; d++) {
    vxi[d] = p.pxi[d] * root;
  }
}

// ----------------------------------------------------------------------
// push_pxi_dt
//
// advance moments according to EM fields

__device__ static void
push_pxi_dt(struct d_particle *p,
	    real exq, real eyq, real ezq, real hxq, real hyq, real hzq)
{
  real dq = p->qni_div_mni * d_dqs;
  real pxm = p->pxi[0] + dq*exq;
  real pym = p->pxi[1] + dq*eyq;
  real pzm = p->pxi[2] + dq*ezq;
  
  real root = dq * rsqrtr(real(1.) + sqr(pxm) + sqr(pym) + sqr(pzm));
  real taux = hxq * root, tauy = hyq * root, tauz = hzq * root;
  
  real tau = real(1.) / (real(1.) + sqr(taux) + sqr(tauy) + sqr(tauz));
  real pxp = ( (real(1.) + sqr(taux) - sqr(tauy) - sqr(tauz)) * pxm
	       +(real(2.)*taux*tauy + real(2.)*tauz)*pym
	       +(real(2.)*taux*tauz - real(2.)*tauy)*pzm)*tau;
  real pyp = ( (real(2.)*taux*tauy - real(2.)*tauz)*pxm
	       +(real(1.) - sqr(taux) + sqr(tauy) - sqr(tauz)) * pym
	       +(real(2.)*tauy*tauz + real(2.)*taux)*pzm)*tau;
  real pzp = ( (real(2.)*taux*tauz + real(2.)*tauy)*pxm
	       +(real(2.)*tauy*tauz - real(2.)*taux)*pym
	       +(real(1.) - sqr(taux) - sqr(tauy) + sqr(tauz))*pzm)*tau;
  
  p->pxi[0] = pxp + dq * exq;
  p->pxi[1] = pyp + dq * eyq;
  p->pxi[2] = pzp + dq * ezq;
}

#define OFF(g, d) o##g[d]
  
#if DIM == DIM_Z

#define INTERPOLATE_FIELD(exq, fldnr, g1, g2)				\
  do {									\
    int ddz = l##g2[2]-l0[2];						\
    exq =								\
      ip_to_grid_m(OFF(g2, 2)) * F3C(fldnr, 0, 0, ddz-1) +		\
      ip_to_grid_0(OFF(g2, 2)) * F3C(fldnr, 0, 0, ddz+0) +		\
      ip_to_grid_p(OFF(g2, 2)) * F3C(fldnr, 0, 0, ddz+1);		\
  } while(0)

#elif DIM == DIM_YZ

#define INTERPOLATE_FIELD(exq, fldnr, g1, g2)				\
  do {									\
    int ddy = l##g1[1]-l0[1], ddz = l##g2[2]-l0[2];			\
    /* printf("C %g [%d,%d,%d]\n", F3C(fldnr, 0, ddy, ddz), 0, ddy, ddz); */ \
    exq =								\
      ip_to_grid_m(OFF(g1, 1)) * ip_to_grid_m(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy-1, ddz-1) +					\
      ip_to_grid_0(OFF(g1, 1)) * ip_to_grid_m(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy+0, ddz-1) +					\
      ip_to_grid_p(OFF(g1, 1)) * ip_to_grid_m(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy+1, ddz-1) +					\
      ip_to_grid_m(OFF(g1, 1)) * ip_to_grid_0(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy-1, ddz+0) +					\
      ip_to_grid_0(OFF(g1, 1)) * ip_to_grid_0(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy+0, ddz+0) +					\
      ip_to_grid_p(OFF(g1, 1)) * ip_to_grid_0(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy+1, ddz+0) +					\
      ip_to_grid_m(OFF(g1, 1)) * ip_to_grid_p(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy-1, ddz+1) +					\
      ip_to_grid_0(OFF(g1, 1)) * ip_to_grid_p(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy+0, ddz+1) +					\
      ip_to_grid_p(OFF(g1, 1)) * ip_to_grid_p(OFF(g2, 2)) *		\
      F3C(fldnr, 0, ddy+1, ddz+1);					\
  } while(0)

#endif

