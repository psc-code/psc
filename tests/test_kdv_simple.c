
#include <mrc_ts.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_params.h>
#include <mrc_bits.h>

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

// FIXME BND
#define CRDX(ix) (MRC_CRDX(crds, (ix)+BND))

static void
fill_ghosts(struct mrc_fld *x, int m_x)
{
  int mx = mrc_fld_dims(x)[0];
  MRC_F1(x, m_x , -2  ) = MRC_F1(x, m_x , mx-2);
  MRC_F1(x, m_x , -1  ) = MRC_F1(x, m_x , mx-1);
  MRC_F1(x, m_x , mx  ) = MRC_F1(x, m_x , 0);
  MRC_F1(x, m_x , mx+1) = MRC_F1(x, m_x , 1);
}

#define Dx(x, m_x, ix)							\
  ((MRC_F1(x, m_x, ix+1) - MRC_F1(x, m_x, ix-1)) / (CRDX(ix+1) - CRDX(ix-1)))

// assumes uniform coordinates!
#define Dxxx(x, m_x, ix)						\
  ((MRC_F1(x, m_x, ix+2) - 2.*MRC_F1(x, m_x, ix+1) + 2.*MRC_F1(x, m_x, ix-1) - MRC_F1(x, m_x, ix-2)) / (2.*powf(CRDX(ix+1) - CRDX(ix), 3.)))

static void
calc_rhs(void *ctx, struct mrc_obj *_rhs, float time, struct mrc_obj *_x)
{
  struct mrc_domain *domain = ctx;
  struct mrc_fld *rhs = (struct mrc_fld *) _rhs, *x = (struct mrc_fld *) _x;
  struct mrc_crds *crds = mrc_domain_get_crds(domain);

  fill_ghosts(x, U);

  mrc_f1_foreach(x, ix, 0, 0) {
    MRC_F1(rhs, U, ix) =  - Dxxx(x, U, ix) + 6. * MRC_F1(x, U, ix) * Dx(x, U, ix);
  } mrc_f1_foreach_end;
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
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_param_int(crds, "sw", BND);
  mrc_crds_set_param_double3(crds, "l", (double[3]) { -8., 0., 0. });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  8., 0., 0. });

  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);

  struct mrc_fld *x = mrc_domain_f1_create(domain);
  mrc_fld_set_name(x, "x");
  mrc_fld_set_param_int(x, "nr_ghosts", BND);
  mrc_fld_set_param_int(x, "nr_comps", NR_FLDS);
  mrc_fld_setup(x);
  mrc_fld_set_comp_name(x, U, "u");

  // setup initial equilibrium and perturbation
  mrc_f1_foreach(x, ix, 0, 0) {
    //    MRC_F1(x, U, ix) = sin(2.*M_PI * CRDX(ix));
    MRC_F1(x, U, ix) = -12. * 1./sqr(cosh(CRDX(ix))); // 3 solitons
  } mrc_f1_foreach_end;

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
