//! \file patchmanager.h This class provides a manager of patches

#include <bitfield3d.h>

struct psc_domainwindow;
struct windowlist_item;

struct psc_patchmanager
{
  struct bitfield3d p1, p2;
  struct bitfield3d *activepatches, *lastactivepatches;
  struct windowlist_item* windows;
  struct psc_case* currentcase;
};

void psc_patchmanager_set_from_options(struct psc_patchmanager *this);

///This creates a patchmanager assuming \a np patches along each axis
///\note Make sure, \a np equals the \a np given to mrc_domain
void psc_patchmanager_setup(struct psc_patchmanager *this, int np[3]);

///Use this to deallocate the patchmanager
void psc_patchmanager_destroy(struct psc_patchmanager *this);

///This function needs to be called every timestep to advance all windows
///and reset the domain if neccessary
void psc_patchmanager_timestep(struct psc_patchmanager* this);

///\brief Creates a new domain-window
///
///Use this function during your case::setup() or afterwards to create new domain-windows
///\param type The type of window to create. Currently, only "movingwindow_z" is supported
struct psc_domainwindow* psc_patchmanager_create_window(struct psc_patchmanager *this, const char* type);