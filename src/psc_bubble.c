
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_bnd_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>

#include <mrc_params.h>

#include <math.h>

struct psc_bubble {
  double BB;
  double nnb;
  double TTe;
  double TTi;
  double MMi;
  double MMach;
  double LLn;
  double LLB;
};

#define to_psc_bubble(psc) mrc_to_subobj(psc, struct psc_bubble)

#define VAR(x) (void *)offsetof(struct psc_bubble, x)
static struct param psc_bubble_descr[] = {
  { "BB"            , VAR(BB)              , PARAM_DOUBLE(.07)    },
  { "nnb"           , VAR(nnb)             , PARAM_DOUBLE(.1)     },
  { "MMi"           , VAR(MMi)             , PARAM_DOUBLE(100.)   },
  { "MMach"         , VAR(MMach)           , PARAM_DOUBLE(3.)     },
  { "LLn"           , VAR(LLn)             , PARAM_DOUBLE(200.)   },
  { "LLB"           , VAR(LLB)             , PARAM_DOUBLE(200./6.)},
  { "TTe"           , VAR(TTe)             , PARAM_DOUBLE(.02)    },
  { "TTi"           , VAR(TTi)             , PARAM_DOUBLE(.02)    },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_bubble_create

static void
psc_bubble_create(struct psc *psc)
{
  // new defaults (dimensionless) for this case
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.nmax = 32000;
  psc->prm.cpum = 5*24.0*60*60;
  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;

  psc->prm.nicell = 50;

  psc->domain.gdims[0] = 2400;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 1600;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part[2] = BND_PART_PERIODIC;

  struct psc_bnd_fields *bnd_fields = 
    psc_push_fields_get_bnd_fields(psc->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "none");

  psc_sort_set_type(psc->sort, "countsort2");
}

// ----------------------------------------------------------------------
// psc_bubble_set_from_options

static void
psc_bubble_set_from_options(struct psc *psc)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);

  psc->domain.length[0] = 3.  * bubble->LLn;
  psc->domain.length[1] = 10. * bubble->LLn; // no y dependence 
  psc->domain.length[2] = 2.  * bubble->LLn;

  psc->domain.corner[0] = -1.5 * bubble->LLn;
  psc->domain.corner[1] = -5.  * bubble->LLn;
  psc->domain.corner[2] = -1.  * bubble->LLn;
}

// ----------------------------------------------------------------------
// psc_bubble_setup_fields

static void
psc_bubble_setup_fields(struct psc *psc, mfields_base_t *flds)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);

  double BB = bubble->BB;
  double LLn = bubble->LLn;
  double LLB = bubble->LLB;
  double MMi = bubble->MMi;
  double MMach = bubble->MMach;
  double TTe = bubble->TTe;

  // FIXME, do we need the ghost points?
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    psc_foreach_3d_g(psc, p, jx, jy, jz) {
      double dx = psc->dx[0], dz = psc->dx[2];
      double xx = CRDX(p, jx), zz = CRDZ(p, jz);

      double x1 = xx;
      double z1 = zz + 0.5*dz + LLn;
      double r1 = sqrt(sqr(x1) + sqr(z1));
      double x2 = xx;
      double z2 = zz + 0.5*dz - LLn;
      double r2 = sqrt(sqr(x2) + sqr(z2));

      if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	F3_BASE(pf, HX, jx,jy,jz) += - BB * sin(M_PI * (LLn - r1)/(2.*LLB)) * z1 / r1;
      }
      if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	F3_BASE(pf, HX, jx,jy,jz) += - BB * sin(M_PI * (LLn - r2)/(2.*LLB)) * z2 / r2;
      }

      x1 = xx + 0.5*dx;
      z1 = zz + LLn;
      r1 = sqrt(sqr(x1) + sqr(z1));
      x2 = xx + 0.5*dx;
      z2 = zz - LLn;
      r2 = sqrt(sqr(x2) + sqr(z2));

      if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	F3_BASE(pf, HZ, jx,jy,jz) += BB * sin(M_PI * (LLn - r1)/(2.*LLB)) * x1 / r1;
      }
      if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	F3_BASE(pf, HZ, jx,jy,jz) += BB * sin(M_PI * (LLn - r2)/(2.*LLB)) * x2 / r2;
      }

      x1 = xx;
      z1 = zz + LLn;
      r1 = sqrt(sqr(x1) + sqr(z1));
      x2 = xx;
      z2 = zz - LLn;
      r2 = sqrt(sqr(x2) + sqr(z2));

      if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	F3_BASE(pf, EY, jx,jy,jz) += MMach * sqrt(TTe/MMi) * BB *
	  sin(M_PI * (LLn - r1)/(2.*LLB)) *
	  sin(M_PI * r1 / LLn);
	F3_BASE(pf, JYI, jx,jy,jz) += BB * M_PI/(2.*LLB) * cos(M_PI * (LLn - r1)/(2.*LLB));
      }
      if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	F3_BASE(pf, EY, jx,jy,jz) += MMach * sqrt(TTe/MMi) * BB *
	  sin(M_PI * (LLn - r2)/(2.*LLB)) *
	  sin(M_PI * r2 / LLn);
	F3_BASE(pf, JYI, jx,jy,jz) += BB * M_PI/(2.*LLB) * cos(M_PI * (LLn - r2)/(2.*LLB));
      }
    } foreach_3d_g_end;
  }
}

// ----------------------------------------------------------------------
// psc_bubble_init_npt

static void
psc_bubble_init_npt(struct psc *psc, int kind, double x[3],
		    struct psc_particle_npt *npt)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);

  double BB = bubble->BB;
  double LLn = bubble->LLn;
  double LLB = bubble->LLB;
  double V0 = bubble->MMach * sqrt(bubble->TTe / bubble->MMi);

  double nnb = bubble->nnb;
  double TTi = bubble->TTi;
  double TTe = bubble->TTe;

  double r1 = sqrt(sqr(x[0]) + sqr(x[2] + LLn));
  double r2 = sqrt(sqr(x[0]) + sqr(x[2] - LLn));

  npt->n = nnb;
  if (r1 < LLn) {
    npt->n += (1. - nnb) * sqr(cos(M_PI / 2. * r1 / LLn));
    if (r1 > 0.0) {
      npt->p[0] += V0 * sin(M_PI * r1 / LLn) * x[0] / r1;
      npt->p[2] += V0 * sin(M_PI * r1 / LLn) * (x[2] + 1.*LLn) / r1;
    }
  }
  if (r2 < LLn) {
    npt->n += (1. - nnb) * sqr(cos(M_PI / 2. * r2 / LLn));
    if (r2 > 0.0) {
      npt->p[0] += V0 * sin(M_PI * r2 / LLn) * x[0] / r2;
      npt->p[2] += V0 * sin(M_PI * r2 / LLn) * (x[2] - 1.*LLn) / r2;
    }
  }

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;

    // electron drift consistent with initial current
    if ((r1 <= LLn) && (r1 >= LLn - 2.*LLB)) {
      npt->p[1] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r1)/(2.*LLB)) / npt->n;
    }
    if ((r2 <= LLn) && (r2 >= LLn - 2.*LLB)) {
      npt->p[1] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r2)/(2.*LLB)) / npt->n;
    }

    npt->T[0] = TTe;
    npt->T[1] = TTe;
    npt->T[2] = TTe;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = bubble->MMi;

    npt->T[0] = TTi;
    npt->T[1] = TTi;
    npt->T[2] = TTi;
    break;
  default:
    assert(0);
  }
}

// ======================================================================
// psc_bubble_ops

struct psc_ops psc_bubble_ops = {
  .name             = "bubble",
  .size             = sizeof(struct psc_bubble),
  .param_descr      = psc_bubble_descr,
  .create           = psc_bubble_create,
  .set_from_options = psc_bubble_set_from_options,
  .setup_fields     = psc_bubble_setup_fields,
  .init_npt         = psc_bubble_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_bubble_ops);

  // psc_create() will create the psc object, create the sub-objects
  // (particle, field pusher and many others) and set the parameter defaults.
  // It will also set the psc_bubble defaults and call psc_bubble_create(),
  // which will change some of the general defaults to match this case.
  struct psc *psc = psc_create(MPI_COMM_WORLD);

  // psc_set_from_options() will override general and bubble psc parameters
  // if given on the command line. It will also call
  // psc_bubble_set_from_options()
  psc_set_from_options(psc);

  // psc_setup() will set up the various sub-objects (particle pusher, ...)
  // and set up the initial domain partition, the particles and the fields.
  // The standard implementation, used here, will set particles using
  // psc_bubble_init_npt and the fields using setup_field()
  psc_setup(psc);

  // psc_view() will just print a whole lot of info about the psc object and
  // sub-objects, in particular all the parameters.
  psc_view(psc);

  // psc_integrate() uses the standard implementation, which does the regular
  // classic PIC time integration loop
  psc_integrate(psc);

  // psc_destroy() just cleans everything up when we're done.
  psc_destroy(psc);

  libmrc_params_finalize();
  MPI_Finalize();
}
