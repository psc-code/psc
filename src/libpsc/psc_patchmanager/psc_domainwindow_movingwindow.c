#include <stdio.h>
#include <math.h>
#include <psc.h>
#include "psc_domainwindow_private.h"
#include <mrc_params.h>

struct psc_movingwindow_z
{
  double speed;		///<speed of the window (physical coordinates)
  double length;	///<length of the window (physical coordinates)
  double position;	///<current position of the lower corner (internal coordinates)
};

#define VAR(x) (void *)offsetof(struct psc_movingwindow_z, x)
static struct param psc_movingwindow_z_descr[] = {
  { "speed"	, VAR(speed)	, PARAM_DOUBLE(1.e-6) },
  { "length"	, VAR(length)	, PARAM_DOUBLE(2.e-6) },
  {}
};
#undef VAR

//Activates a whole xy-subdomain at position [z = slice]
static void set_slice(struct bitfield3d* patches, int slice, int value)
{
  printf("set slice %d to %d\n", slice, value);
  for(int y=0; y<ppsc->domain.np[1]; ++y)
  {
    for(int x=0; x<ppsc->domain.np[0]; ++x)
    {
      bitfield3d_set(patches, (int[3]){x,y,slice}, value);
    }
  }
}

static void psc_movingwindow_z_setup(struct psc_domainwindow* this)
{
  struct psc_movingwindow_z* window = mrc_to_subobj(this, struct psc_movingwindow_z);
  window->position = 0.;
  
  //Activate all patches that are touched by the initial window
  double patchlength = ppsc->domain.length[2] / ppsc->domain.np[2];	//The z-length of one patch in physical coordinates
  double windowlength = window->length;		//The z-length of the window in physical coordinates
  int npatches = ceil(windowlength / patchlength);

  for(int i=0; i<npatches; ++i)
  {
    if(i < ppsc->domain.np[2]) set_slice(&this->activepatches, i, 1);
  }
}

inline int sgn(int x)
{
  return x < 0 ? -1 : 1;
}

static void psc_movingwindow_z_timestep(struct psc_domainwindow* this, int timestep, double t)
{
  struct psc_movingwindow_z* window = mrc_to_subobj(this, struct psc_movingwindow_z);
  
  double patchlength = ppsc->domain.length[2] / ppsc->domain.np[2];	//The z-length of one patch in physical coordinates
  double v = window->speed / ppsc->coeff.wl;			//velocity in physical coordinates
  
  //Deactivate obsolete slices
  //lowest active slice
  int ls = window->position / patchlength;
  int ls_new = (v * t) / patchlength;
  
  while(ls_new != ls)
  {
    if(ls < ppsc->domain.np[2]) set_slice(&this->activepatches, ls, 0);
    ls += sgn(ls_new - ls);
  }
  
  //Activate new slices
  double windowlength = window->length;			//The z-length of the window in physical coordinates
  int hs = ceil((window->position + windowlength) / patchlength);
  int hs_new = ceil((v * t + windowlength) / patchlength);
  
  while(hs_new != hs)
  {
    if(hs < ppsc->domain.np[2]) set_slice(&this->activepatches, hs, 1);
    hs += sgn(hs_new - hs);
  }
  
  window->position = v * t;
}

struct psc_domainwindow_ops psc_movingwindow_z_ops = {
  .name             = "movingwindow_z",
  .size             = sizeof(struct psc_movingwindow_z),
  .setup            = psc_movingwindow_z_setup,
  .param_descr      = psc_movingwindow_z_descr,
  .timestep         = psc_movingwindow_z_timestep,
};
