
#include "mrc_ggcm_gridx_gen.h"

#include <mrc_domain.h>
#include <mrc_params.h>
#include <math.h>
#include <assert.h>

#define ITERMAX 2000

// ----------------------------------------------------------------------
// integrate_rk4

static double
integrate_rk4(struct mrc_crds_gen *gen, double xi, double dxi, double fak,
              double (*dx_func)(struct mrc_crds_gen *gen, double x, double fak))
{
  int ns = fabs((xi - 1.) / dxi + .5) + 5;
  double ddxi = (xi - 1.) / ns;

  double x = gen->xl;
  for (int i = 0; i < ns; i++) {
    // 4th order runge kutta
    double xk1 = ddxi * dx_func(gen, x, fak);
    double xk2 = ddxi * dx_func(gen, x + 0.5 * xk1, fak);
    double xk3 = ddxi * dx_func(gen, x + 0.5 * xk2, fak);
    double xk4 = ddxi * dx_func(gen, x + xk3, fak);
    x += (1.0/6.0) * (xk1 + xk4) + (1.0/3.0) * (xk2 + xk3);
  }

  return x;
}

// ----------------------------------------------------------------------
// generate_ggcm_x_grid

void
generate_ggcm_x_grid(struct mrc_crds_gen *gen, double *xx, double *dx,
                     double (*dx_func)(struct mrc_crds_gen *gen, double x, double fak))
{
  double fak, dxi, s;

  //  printf("gridx: n = %d\n", gen->n);

  // convergence test xni, the number of steps in the rk4 integrator
  double xn1 = 0.0;
  double xn1o = 0.0;
  double rel = 0.0;

  dxi = 0.1;
  fak = 1.0;
  for (int k = 0; k < ITERMAX; k++) {
    xn1 = integrate_rk4(gen, gen->n, dxi, fak, dx_func);
    rel = fabs(xn1 - xn1o) / (fabs(xn1) + fabs(xn1o));
    // printf("%d: convergence test: %.18g, %.18g, %.18g, %.18g\n", k, dxi, xn1, xn1o, rel);
    if (rel < 1e-8)
      break;
    xn1o = xn1;
    dxi *= 0.95;
  }
  //  printf("gridx: convergence test: %g, %g\n", dxi, xn1);

  // find fak parameter that gives a grid with the correct length
  fak = 0.1;
  s = 0.1;
  for (int k = 0; k < ITERMAX; k++) {
    xn1 = integrate_rk4(gen, gen->n, dxi, fak, dx_func);
    //    printf("%d: fak=%g xn1=%g xn=%g\n", k, fak, xn1, gen->xh);

    if (xn1 > gen->xh) {
      if (fabs(xn1 - gen->xh) < 1e-6)
	      break;
      fak -= s;
      s *= 0.1;
    } else {
      fak += s;
    }
  }
  //  printf("gridx: fak=%g xn1=%g xn=%g\n", fak, xn1, gen->xh);

  double length, target_length;
  length = xn1 - gen->xl;
  target_length = gen->xh - gen->xl;
  if (fabs(length - target_length) / target_length > 0.01) {
    mpi_printf(MPI_COMM_WORLD,
               "Error: grid did not converge to within 1%% of target length!\n");
    mpi_printf(MPI_COMM_WORLD, "(length = %g, target = %g)\n", length, target_length);
    assert(0);
  }

  int step_out = gen->n / 20;
  if (step_out > 20) step_out = 20;

  for (int i = -gen->sw; i < gen->n + gen->sw; i++) {
    xx[i] = integrate_rk4(gen, i+1, dxi, fak, dx_func);
    dx[i] = dx_func(gen, xx[i], fak);

    // if (i % step_out == 1) {
    //   printf("%d: xx=%g dx=%g\n", i, xx[i], dx[i]);
    // }
  }
}
