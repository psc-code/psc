
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_collision.h>
#include <psc_balance.h>

#include <psc_fields_as_c.h>
#include <psc_fields_single.h>
#include "fields.hxx"

#include <mrc_params.h>

#include <math.h>

using Fields = Fields3d<fields_t>;

struct psc_test_em_wave {
  // params
  double k[3];
  double B[3];
  double tol;

  // state
  double om;
  double E[3];
};

#define psc_test_em_wave(psc) mrc_to_subobj(psc, struct psc_test_em_wave)

#define VAR(x) (void *)offsetof(struct psc_test_em_wave, x)
static struct param psc_test_em_wave_descr[] = {
  { "k"               , VAR(k)                 , PARAM_DOUBLE3(2., 0., 1.)  },
  { "B"               , VAR(B)                 , PARAM_DOUBLE3(0., 1., 0.)  },
  { "tol"             , VAR(tol)               , PARAM_DOUBLE(1e-6)         },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_test_em_wave_create

static void
psc_test_em_wave_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nicell = 1;
  psc->prm.cfl = 1.;

  psc->domain.gdims[0] = 16;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 16;

  psc->domain.length[0] = 2. * M_PI;
  psc->domain.length[1] = 2. * M_PI;
  psc->domain.length[2] = 2. * M_PI;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
 
  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;

  psc_push_particles_set_type(psc->push_particles, "1vbec_double");
}

// ----------------------------------------------------------------------
// psc_test_em_wave_setup

static void
psc_test_em_wave_setup(struct psc *psc)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);

  double *k = sub->k, *B = sub->B, *E = sub->E;

  assert(k[0]*B[0] + k[1]*B[1] + k[2]*B[2] == 0);
  double k_norm = sqrt(sqr(k[0]) + sqr(k[1]) + sqr(k[2]));
  sub->om = /* c * */ k_norm;
  E[0] = 1./sub->om * (k[1] * B[2] - k[2] * B[1]);
  E[1] = 1./sub->om * (k[2] * B[0] - k[0] * B[2]);
  E[2] = 1./sub->om * (k[0] * B[1] - k[1] * B[0]);

  psc_setup_super(psc);

}

// ----------------------------------------------------------------------
// psc_test_em_wave_init_field

static double
psc_test_em_wave_init_field(struct psc *psc, double crd[3], int m)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);
  double *k = sub->k, *B = sub->B, *E = sub->E;
  double x = crd[0], y = crd[1], z = crd[2];

  switch (m) {
  case EX: return E[0] * sin(-(k[0]*x + k[1]*y + k[2]*z));
  case EY: return E[1] * sin(-(k[0]*x + k[1]*y + k[2]*z));
  case EZ: return E[2] * sin(-(k[0]*x + k[1]*y + k[2]*z));
  case HX: return B[0] * sin(-(k[0]*x + k[1]*y + k[2]*z));
  case HY: return B[1] * sin(-(k[0]*x + k[1]*y + k[2]*z));
  case HZ: return B[2] * sin(-(k[0]*x + k[1]*y + k[2]*z));
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// is_equal_tol

static int
is_equal_tol(fields_t::real_t tol, fields_t::real_t v1, fields_t::real_t v2)
{
  if (fabs(v1 - v2) < tol) {
    return 0;
  }
  printf("is_equal_tol: v1 %g v2 %g diff %g (tol %g)\n", v1, v2, v1-v2, tol);
  return 1;
}

// ----------------------------------------------------------------------
// psc_test_em_wave_check

static void
psc_test_em_wave_check(struct psc *psc)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);
  double om = -sub->om; // FIXME, looks like we need a minus sign here, but why?
  double t = psc->timestep * psc->dt;
  double *k = sub->k, *E = sub->E, *B = sub->B;
  double tol = sub->tol;

  printf("=== checking EM fields at time step %d\n", psc->timestep);
  mfields_t mf = psc->flds->get_as<mfields_t>(EX, EX + 6);

  int failed = 0;
  for (int p = 0; p < mf.n_patches(); p++) {
    Fields F(mf[p]);

    foreach_3d(psc, p, jx,jy,jz, 0, 0) {
      double dx = psc->patch[p].dx[0], dy = psc->patch[p].dx[1], dz = psc->patch[p].dx[2];
      double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);

      failed += is_equal_tol(tol, F(EX, jx,jy,jz),
			     E[0] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy         ) + k[2]*(zz         ))));
      failed += is_equal_tol(tol, F(EY, jx,jy,jz),
			     E[1] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy + .5f*dy) + k[2]*(zz         ))));
      failed += is_equal_tol(tol, F(EZ, jx,jy,jz),
			     E[2] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy         ) + k[2]*(zz + .5f*dz))));

      failed += is_equal_tol(tol, F(HX, jx,jy,jz),
			     B[0] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy + .5f*dy) + k[2]*(zz + .5f*dz))));
      failed += is_equal_tol(tol, F(HY, jx,jy,jz),
			     B[1] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy         ) + k[2]*(zz + .5f*dz))));
      failed += is_equal_tol(tol, F(HZ, jx,jy,jz),
			     B[2] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy + .5f*dy) + k[2]*(zz         ))));
    } foreach_3d_end;
  }
  assert(failed == 0);

  mf.put_as(psc->flds, 0, 0);
}

// ----------------------------------------------------------------------
// psc_test_em_wave_check_single

static void
psc_test_em_wave_check_single(struct psc *psc)
{
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);
  double om = -sub->om; // FIXME, looks like we need a minus sign here, but why?
  double t = psc->timestep * psc->dt;
  double *k = sub->k, *E = sub->E, *B = sub->B;
  double tol = sub->tol;

  printf("=== single: checking EM fields at time step %d\n", psc->timestep);
  mfields_single_t mf = psc->flds->get_as<mfields_single_t>(EX, EX + 6);

  int failed = 0;
  for (int p = 0; p < mf.n_patches(); p++) {
    Fields3d<fields_single_t> F(mf[p]);

    foreach_3d(psc, p, jx,jy,jz, 0, 0) {
      double dx = psc->patch[p].dx[0], dy = psc->patch[p].dx[1], dz = psc->patch[p].dx[2];
      double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);

      failed += is_equal_tol(tol, F(EX, jx,jy,jz),
			     E[0] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy         ) + k[2]*(zz         ))));
      failed += is_equal_tol(tol, F(EY, jx,jy,jz),
			     E[1] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy + .5f*dy) + k[2]*(zz         ))));
      failed += is_equal_tol(tol, F(EZ, jx,jy,jz),
			     E[2] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy         ) + k[2]*(zz + .5f*dz))));

      failed += is_equal_tol(tol, F(HX, jx,jy,jz),
			     B[0] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy + .5f*dy) + k[2]*(zz + .5f*dz))));
      failed += is_equal_tol(tol, F(HY, jx,jy,jz),
			     B[1] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy         ) + k[2]*(zz + .5f*dz))));
      failed += is_equal_tol(tol, F(HZ, jx,jy,jz),
			     B[2] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy + .5f*dy) + k[2]*(zz         ))));
    } foreach_3d_end;
  }
  assert(failed == 0);

  mf.put_as(psc->flds, 0, 0);
}

// ----------------------------------------------------------------------
// psc_test_em_wave_check_vpic

static void
psc_test_em_wave_check_vpic(struct psc *psc)
{
#if 0
  struct psc_test_em_wave *sub = psc_test_em_wave(psc);
  double om = -sub->om; // FIXME, looks like we need a minus sign here, but why?
  double t = psc->timestep * psc->dt;
  double *k = sub->k, *E = sub->E, *B = sub->B;
  double tol = sub->tol;
#endif

  printf("=== vpic: checking EM fields at time step %d\n", psc->timestep);
  mfields_single_t mf = psc->flds->get_as<mfields_single_t>(EX, EX + 6);
  
#if 0
  int failed = 0;
  for (int p = 0; p < mf.nr_patches(); p++) {
    fields_single_t flds = mfields_single_t(mflds)[p];

    foreach_3d(psc, p, jx,jy,jz, 0, 0) {
      double dx = psc->patch[p].dx[0], dy = psc->patch[p].dx[1], dz = psc->patch[p].dx[2];
      double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);

      failed += is_equal_tol(tol, F(EX, jx,jy,jz),
			     E[0] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy         ) + k[2]*(zz         ))));
      failed += is_equal_tol(tol, F(EY, jx,jy,jz),
			     E[1] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy + .5f*dy) + k[2]*(zz         ))));
      failed += is_equal_tol(tol, F(EZ, jx,jy,jz),
			     E[2] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy         ) + k[2]*(zz + .5f*dz))));

      failed += is_equal_tol(tol, F(HX, jx,jy,jz),
			     B[0] * sin(om*t - (k[0]*(xx         ) + k[1]*(yy + .5f*dy) + k[2]*(zz + .5f*dz))));
      failed += is_equal_tol(tol, F(HY, jx,jy,jz),
			     B[1] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy         ) + k[2]*(zz + .5f*dz))));
      failed += is_equal_tol(tol, F(HZ, jx,jy,jz),
			     B[2] * sin(om*t - (k[0]*(xx + .5f*dx) + k[1]*(yy + .5f*dy) + k[2]*(zz         ))));
    } foreach_3d_end;
  }
  assert(failed == 0);
#endif

  mf.put_as(psc->flds, 0, 0);
}

// ----------------------------------------------------------------------
// psc_test_em_wave_step

static void
psc_test_em_wave_step(struct psc *psc)
{
  psc_output(psc);

  psc_test_em_wave_check(psc);
  psc_test_em_wave_check_single(psc);
  psc_test_em_wave_check_vpic(psc);
  
  psc_push_fields_push_H(psc->push_fields, psc->flds, .5);
  psc_push_fields_step_b2(psc->push_fields, psc->flds);
  psc_push_fields_step_a(psc->push_fields, psc->flds);
}

// ======================================================================
// psc_test_em_wave_ops

struct psc_ops psc_test_em_wave_ops = {
  .name             = "test_em_wave",
  .size             = sizeof(struct psc_test_em_wave),
  .param_descr      = psc_test_em_wave_descr,
  .create           = psc_test_em_wave_create,
  .setup            = psc_test_em_wave_setup,
  .init_field       = psc_test_em_wave_init_field,
  .step             = psc_test_em_wave_step,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_em_wave_ops);
}
