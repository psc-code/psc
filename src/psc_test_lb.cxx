
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_bnd_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>

#include <mrc_params.h>

#include <math.h>

struct psc_test_lb {
  double nb;
  double nn;
  double sy;
  double sz;
  double sigma;
};

#define psc_test_lb(psc) mrc_to_subobj(psc, struct psc_test_lb)

#define VAR(x) (void *)offsetof(struct psc_test_lb, x)
static struct param psc_test_lb_descr[] = {
  { "nb"            , VAR(nb)              , PARAM_DOUBLE(.1)     },
  { "nn"            , VAR(nn)              , PARAM_DOUBLE(1.)     },
  { "sy"            , VAR(sy)              , PARAM_DOUBLE(1.)     },
  { "sz"            , VAR(sz)              , PARAM_DOUBLE(1.)     },
  { "sigma"         , VAR(sigma)           , PARAM_DOUBLE(1.)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_test_lb_create

static void
psc_test_lb_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 21;
  psc->prm.nicell = 50;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 128;
  psc->domain.gdims[2] = 128;

  psc->domain.length[0] = 1.;
  psc->domain.length[1] = 2.;
  psc->domain.length[2] = 2.;

  psc->domain.corner[0] = 0.;
  psc->domain.corner[1] = -1.;
  psc->domain.corner[2] = -1.;

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

  struct psc_bnd_fields *bnd_fields = 
    psc_push_fields_get_bnd_fields(psc->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "none");
}

// ----------------------------------------------------------------------
// psc_test_lb_init_npt

static void
psc_test_lb_init_npt(struct psc *psc, int kind, double x[3],
		    struct psc_particle_npt *npt)
{
  struct psc_test_lb *sub = psc_test_lb(psc);

  double d = (x[1] * sub->sy + x[2] * sub->sz) / sqrt(sqr(sub->sy) + sqr(sub->sz));
  npt->n = sub->nb + sub->nn * exp(-sqr(d / sub->sigma));
}

// ======================================================================
// psc_test_lb_ops

struct psc_ops psc_test_lb_ops = {
  .name             = "bubble",
  .size             = sizeof(struct psc_test_lb),
  .param_descr      = psc_test_lb_descr,
  .create           = psc_test_lb_create,
  .init_npt         = psc_test_lb_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_lb_ops);
}
