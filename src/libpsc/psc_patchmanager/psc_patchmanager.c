#include <mrc_domain.h>
#include <mrc_domain_dynamic.h>

#include "psc.h"
#include "psc_balance.h"

#include "psc_bnd_fields.h"
#include "psc_push_fields.h"
#include "psc_bnd.h"

#include "psc_domainwindow_private.h"

struct windowlist_item
{
  struct psc_domainwindow* v;
  struct windowlist_item* next;
};

/* TODO
 * Add non-existing patches to the load-calculation in psc_balance
 * As of now, they'll only be counted at the next LB-step
 * */

void psc_patchmanager_set_from_options(struct psc_patchmanager* this)
{
  unsigned int* np = (unsigned*) ppsc->domain.np;
  bitfield3d_create(&this->p1, np);
  bitfield3d_create(&this->p2, np);
  
  this->activepatches = &this->p1;
  this->lastactivepatches = &this->p2;
  
  this->windows = NULL;
  this->currentcase = NULL;	//TODO: Where should I set this elegantly?
}

void psc_patchmanager_setup(struct psc_patchmanager *this, int np[3])
{
  //Setup all windows
  for(struct windowlist_item* iter = this->windows; iter != NULL; iter = iter->next)
  {
    psc_domainwindow_setup(iter->v);
  }
}

void psc_patchmanager_destroy(struct psc_patchmanager *this)
{
  bitfield3d_destroy(&this->p1);
  bitfield3d_destroy(&this->p2);
}

void psc_patchmanager_timestep(struct psc_patchmanager* this)
{
  //Swap activepatches with lastactivepatches
  struct bitfield3d* b = this->lastactivepatches;
  this->lastactivepatches = this->activepatches;
  this->activepatches = b;
  //Clear list of active patches
  bitfield3d_fill(this->activepatches, 0);
  //Iterate over all windows and update the activepatch-list
  double t = ppsc->dt * ppsc->timestep;
  
  for(struct windowlist_item* iter = this->windows; iter != NULL; iter = iter->next)
  {
    psc_domainwindow_timestep(iter->v, ppsc->timestep, t);
    bitfield3d_merge(this->activepatches, &iter->v->activepatches);
  }
  
  //Invalidate load-balancer (only if the domain changes)
  if( bitfield3d_compare(this->activepatches, this->lastactivepatches) != 0 )
  {
    psc_balance_set_param_int(ppsc->balance, "force_update", 1);
  }
}

struct psc_domainwindow* psc_patchmanager_create_window(struct psc_patchmanager *this, const char* type)
{
  struct psc_domainwindow* window = psc_domainwindow_create(MPI_COMM_WORLD);
  psc_domainwindow_set_type(window, type);
  psc_domainwindow_set_param_int3(window, "np", ppsc->domain.np);
  psc_domainwindow_set_from_options(window);
  
  //Add the window to the list
  struct windowlist_item* listitem = malloc(sizeof(struct windowlist_item));
  listitem->v = window;
  listitem->next = this->windows;
  this->windows = listitem;
  
  return window;
}
