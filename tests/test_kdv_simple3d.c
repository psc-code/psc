
#include <mrc_ts.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_params.h>
#include <mrc_bits.h>
#include <mrc_ddc.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================

#define BND 2

enum {
  U,
  NR_FLDS,
};

// ======================================================================

#define CRDX(ix) (MRC_CRDX(crds, (ix)))

#define Dx(x, m_x, ix,iy,iz)						\
  ((MRC_F3(x, m_x, ix+1,iy,iz) - MRC_F3(x, m_x, ix-1,iy,iz))		\
   / (CRDX(ix+1) - CRDX(ix-1)))

// assumes uniform coordinates!
#define Dxxx(x, m_x, ix,iy,iz)						\
  ((MRC_F3(x, m_x, ix+2,iy,iz) - 2.*MRC_F3(x, m_x, ix+1,iy,iz) +	\
    2.*MRC_F3(x, m_x, ix-1,iy,iz) - MRC_F3(x, m_x, ix-2,iy,iz))		\
   / (2.*powf(CRDX(ix+1) - CRDX(ix), 3.)))

static void
calc_rhs(void *ctx, struct mrc_obj *_rhs, float time, struct mrc_obj *_x)
{
  struct mrc_domain *domain = ctx;
  struct mrc_fld *rhs = (struct mrc_fld *) _rhs, *x = (struct mrc_fld *) _x;
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  struct mrc_ddc *ddc = mrc_domain_get_ddc(domain);
  
  mrc_ddc_fill_ghosts_fld(ddc, 0, NR_FLDS, x);

  mrc_fld_foreach(x, ix, iy, iz, 0, 0) {
    MRC_F3(rhs, U, ix,iy,iz) =  - Dxxx(x, U, ix,iy,iz) +
      6. * MRC_F3(x, U, ix,iy,iz) * Dx(x, U, ix,iy,iz);
  } mrc_fld_foreach_end;
}

// ======================================================================

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);

  // set defaults
  mrc_domain_set_param_int3(domain, "m", (int [3]) { 160, 1, 1 });
  mrc_domain_set_param_int(domain, "bcx", BC_PERIODIC);
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_param_int(crds, "sw", BND);
  mrc_crds_set_param_double3(crds, "l", (double[3]) { -8., 0., 0. });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  8., 1., 1. });

  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);

  struct mrc_fld *x = mrc_domain_fld_create(domain, BND, "u");
  mrc_fld_set_name(x, "x");
  mrc_fld_setup(x);

  // setup initial equilibrium and perturbation
  mrc_fld_foreach(x, ix, iy, iz, 0, 0) {
    MRC_F3(x, U, ix,iy,iz) = -12. * 1./sqr(cosh(CRDX(ix))); // 3 solitons
  } mrc_fld_foreach_end;

  // run time integration
  struct mrc_ts *ts = mrc_ts_create_std(MPI_COMM_WORLD, NULL, NULL);
  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(x));
  mrc_ts_set_rhs_function(ts, calc_rhs, domain);
  mrc_ts_set_from_options(ts);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  mrc_ts_view(ts);
  mrc_ts_destroy(ts);

  mrc_fld_destroy(x);

  MPI_Finalize();
  return 0;
}
