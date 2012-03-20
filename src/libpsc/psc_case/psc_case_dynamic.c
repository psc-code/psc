#include "psc.h"		//Required for global parameters etc
#include "psc_case_private.h"	//Interfaces for a case
#include "psc_pulse.h"		//Interfaces for laser-pulses

#include "psc_domainwindow.h"
#include "psc_push_fields.h"
#include "psc_bnd_fields.h"	//Also required for laser-pulses

#include <mrc_params.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>


//This contains all parameters we need
struct psc_case_dynamic {
  // parameters
  double radius;         //< radius of spherical shell
  int nmax;
};

//Expose all parameters to the command-line parser
#define VAR(x) (void *)offsetof(struct psc_case_dynamic, x)
static struct param psc_case_dynamic_descr[] = {
  { "radius"        , VAR(radius)          , PARAM_DOUBLE(2.e-6)              },
  { "nmax", VAR(nmax), PARAM_INT(201) },
  {}
};
#undef VAR

static void
psc_case_dynamic_set_from_options(struct psc_case *_case)
{
  struct psc_case_dynamic *dynamic = mrc_to_subobj(_case, struct psc_case_dynamic);

  //Here we can override PSCs default parameters
  ppsc->prm.nicell = 5;	//number of particles in cell to generate density 1
  ppsc->prm.nmax = dynamic->nmax;	//Number of timesteps

  //Set the domain's physical size
  ppsc->domain.length[0] = 10 * 1e-6;
  ppsc->domain.length[1] = 10 * 1e-6;
  ppsc->domain.length[2] = 100 * 1e-6;

  //Set the grid size
  ppsc->domain.gdims[0] = 1;
  ppsc->domain.gdims[1] = 1;
  ppsc->domain.gdims[2] = 1024;
  
  ppsc->domain.np[0] = 1;
  ppsc->domain.np[1] = 1;
  ppsc->domain.np[2] = 32;
  
  ppsc->use_dynamic_patches = true;	//This will enable us to later use the patchmanager

  //Set boundary conditions for the fields
  ppsc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_lo[2] = BND_FLD_OPEN;
  ppsc->domain.bnd_fld_hi[2] = BND_FLD_OPEN;
  //Set boundary consitions for the particles

  ppsc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;

  //Insert a laser pulse
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse_z1 = psc_bnd_fields_get_pulse_z1(bnd_fields);
  psc_pulse_set_type(pulse_z1, "gauss");
  psc_pulse_set_param_double3(pulse_z1, "m",  (double[3]) { .5 * ppsc->domain.length[0], .5 * ppsc->domain.length[1], 5e-6});
  psc_pulse_set_param_double3(pulse_z1, "dm", (double[3]) { 2e-6, 2e-6, 2e-6 });
  psc_pulse_set_param_double(pulse_z1, "amplitude_s", 1.);
}

static void
psc_case_dynamic_init_npt(struct psc_case *_case, int kind, double x[3],
			      struct psc_particle_npt *npt)
{
  //Downcast psc_case to psc_case_dynamic
  //struct psc_case_dynamic *dynamiccase = mrc_to_subobj(_case, struct psc_case_dynamic);
  
  //Calculate the density
  double xr[3] = {0.};	//position in physical coordinates centered in the computational domain
  for (int d = 0; d < 3; d++) {
    if(ppsc->domain.gdims[d] > 1)
    {
      xr[d] = x[d] * ppsc->coeff.ld;
      xr[d] -= ppsc->domain.length[d] / 2.;
    }
  };
  
  xr[2] += ppsc->domain.length[2] * 0.35;
  double r = sqrt(xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]);
  //double dens = r < (1e-6)*(1e-6) ? 0.1 : 0.0;
  double dens = 0;
  //if(r > 2e-6 ) 
  dens = 0.5 * (1.-fabs(((r) / 2e-6)));
  if(dens < 0.) dens = 0.;

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1;
    npt->n = dens;
    npt->p[2] = 0.;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = 100;
    npt->n = 0;//dens;
    break;
  default:
    assert(0);
  }
}

void psc_case_dynamic_setup(struct psc_case *_case)
{
  struct psc_domainwindow* window = psc_patchmanager_create_window(&ppsc->patchmanager, "movingwindow_z");
  psc_domainwindow_set_param_double(window, "length", 10e-6);
  psc_domainwindow_set_param_double(window, "speed", ppsc->prm.cc );
}

struct psc_case_ops psc_case_dynamic_ops = {
  .name             = "dynamic",
  .size             = sizeof(struct psc_case_dynamic),
  .param_descr      = psc_case_dynamic_descr,
  .set_from_options = psc_case_dynamic_set_from_options,
  .setup            = psc_case_dynamic_setup,
  .init_npt         = psc_case_dynamic_init_npt,
};
