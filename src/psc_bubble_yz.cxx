
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
  double nn0;
  double MMach;
  double LLn;
  double LLB;
  double LLz;
  double LLy;
  double TTe;
  double TTi;
  double MMi;
};

#define to_psc_bubble(psc) mrc_to_subobj(psc, struct psc_bubble)

#define VAR(x) (void *)offsetof(struct psc_bubble, x)
static struct param psc_bubble_descr[] = {
  { "BB"            , VAR(BB)              , PARAM_DOUBLE(.07)    },
  { "nnb"           , VAR(nnb)             , PARAM_DOUBLE(.1)     },
  { "nn0"           , VAR(nn0)             , PARAM_DOUBLE(1.)     },
  { "MMach"         , VAR(MMach)           , PARAM_DOUBLE(3.)     },
  { "LLn"           , VAR(LLn)             , PARAM_DOUBLE(200.)   },
  { "LLB"           , VAR(LLB)             , PARAM_DOUBLE(200./6.)},
  { "LLz"           , VAR(LLz)             , PARAM_DOUBLE(0.)     },
  { "LLy"           , VAR(LLy)             , PARAM_DOUBLE(0.)     },
  { "TTe"           , VAR(TTe)             , PARAM_DOUBLE(.02)    },
  { "TTi"           , VAR(TTe)             , PARAM_DOUBLE(.02)    },
  { "MMi"           , VAR(MMi)             , PARAM_DOUBLE(100.)   },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_bubble_create

static void
psc_bubble_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 32000;
  psc->prm.nicell = 10;

  auto grid_params = GridParams{};
  grid_params.gdims = { 1, 64, 256 };

  grid_params.bnd_fld_lo = { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC };
  grid_params.bnd_fld_hi = { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC };
  grid_params.bnd_part_lo = { BND_PART_PERIODIC, BND_PART_PERIODIC, BND_PART_PERIODIC };
  grid_params.bnd_part_hi = { BND_PART_PERIODIC, BND_PART_PERIODIC, BND_PART_PERIODIC };

  psc->domain_ = grid_params;

  struct psc_bnd_fields *bnd_fields = 
    psc_push_fields_get_bnd_fields(psc->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "none");
}

// ----------------------------------------------------------------------
// psc_bubble_setup

static void
psc_bubble_setup(struct psc *psc)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);
  
  if (bubble->LLy == 0.) {
    bubble->LLy = 2. * bubble->LLn;
  }
  if (bubble->LLz == 0.) {
    bubble->LLz = 3. * bubble->LLn;
  }
  psc->domain_.length = { bubble->LLn, bubble->LLy, bubble->LLz };
  psc->domain_.corner = { 0, -.5 * bubble->LLy, -.5 * bubble->LLz };

  psc_setup_super(psc);

  MPI_Comm comm = psc_comm(psc);
  mpi_printf(comm, "lambda_D = %g\n", sqrt(bubble->TTe));
}

// ----------------------------------------------------------------------
// psc_bubble_read

static void
psc_bubble_read(struct psc *psc, struct mrc_io *io)
{
  psc_read_super(psc, io);
}

// ----------------------------------------------------------------------
// psc_bubble_init_field

static double
psc_bubble_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);

  double BB = bubble->BB;
  double LLn = bubble->LLn;
  double LLy = bubble->LLy;
  double LLB = bubble->LLB;
  double MMi = bubble->MMi;
  double MMach = bubble->MMach;
  double TTe = bubble->TTe;

  double z1 = x[2];
  double y1 = x[1] + .5 * LLy;
  double r1 = sqrt(sqr(z1) + sqr(y1));
  double z2 = x[2];
  double y2 = x[1] - .5 * LLy;
  double r2 = sqrt(sqr(z2) + sqr(y2));

  double rv = 0.;
  switch (m) {
  case HZ:
    if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
      rv += - BB * sin(M_PI * (LLn - r1)/(2.*LLB)) * y1 / r1;
    }
    if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
      rv += - BB * sin(M_PI * (LLn - r2)/(2.*LLB)) * y2 / r2;
    }
    return rv;

  case HY:
    if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
      rv += BB * sin(M_PI * (LLn - r1)/(2.*LLB)) * z1 / r1;
    }
    if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
      rv += BB * sin(M_PI * (LLn - r2)/(2.*LLB)) * z2 / r2;
    }
    return rv;

  case EX:
    if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
      rv += MMach * sqrt(TTe/MMi) * BB *
	sin(M_PI * (LLn - r1)/(2.*LLB)) * sin(M_PI * r1 / LLn);
    }
    if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
      rv += MMach * sqrt(TTe/MMi) * BB *
	sin(M_PI * (LLn - r2)/(2.*LLB)) * sin(M_PI * r2 / LLn);
    }
    return rv;

  case JXI:
    if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
      rv += BB * M_PI/(2.*LLB) * cos(M_PI * (LLn - r1)/(2.*LLB));
    }
    if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
      rv += BB * M_PI/(2.*LLB) * cos(M_PI * (LLn - r2)/(2.*LLB));
    }
    return rv;

  default:
    return 0.;
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
  double LLy = bubble->LLy;
  double LLn = bubble->LLn;
  double LLB = bubble->LLB;
  double V0 = bubble->MMach * sqrt(bubble->TTe / bubble->MMi);

  double nnb = bubble->nnb;
  double nn0 = bubble->nn0;

  double TTe = bubble->TTe, TTi = bubble->TTi;

  double r1 = sqrt(sqr(x[2]) + sqr(x[1] + .5 * LLy));
  double r2 = sqrt(sqr(x[2]) + sqr(x[1] - .5 * LLy));

  npt->n = nnb;
  if (r1 < LLn) {
    npt->n += (nn0 - nnb) * sqr(cos(M_PI / 2. * r1 / LLn));
    if (r1 > 0.0) {
      npt->p[2] += V0 * sin(M_PI * r1 / LLn) * x[2] / r1;
      npt->p[1] += V0 * sin(M_PI * r1 / LLn) * (x[1] + .5 * LLy) / r1;
    }
  }
  if (r2 < LLn) {
    npt->n += (nn0 - nnb) * sqr(cos(M_PI / 2. * r2 / LLn));
    if (r2 > 0.0) {
      npt->p[2] += V0 * sin(M_PI * r2 / LLn) * x[2] / r2;
      npt->p[1] += V0 * sin(M_PI * r2 / LLn) * (x[1] - .5 * LLy) / r2;
    }
  }

  switch (kind) {
  case 0: // electrons
    // electron drift consistent with initial current
    if ((r1 <= LLn) && (r1 >= LLn - 2.*LLB)) {
      npt->p[0] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r1)/(2.*LLB)) / npt->n;
    }
    if ((r2 <= LLn) && (r2 >= LLn - 2.*LLB)) {
      npt->p[0] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r2)/(2.*LLB)) / npt->n;
    }

    npt->T[0] = TTe;
    npt->T[1] = TTe;
    npt->T[2] = TTe;
    break;
  case 1: // ions
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

struct psc_ops_bubble : psc_ops {
  psc_ops_bubble() {
    name             = "bubble";
    size             = sizeof(struct psc_bubble);
    param_descr      = psc_bubble_descr;
    create           = psc_bubble_create;
    setup            = psc_bubble_setup;
    read             = psc_bubble_read;
    init_field       = psc_bubble_init_field;
    init_npt         = psc_bubble_init_npt;
  }
} psc_bubble_ops;

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_bubble_ops);
}
