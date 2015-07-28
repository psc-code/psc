
#include "psc_cuda2.h"

#include "psc_fields_as_single.h"

// ----------------------------------------------------------------------
// cuda2_push_flds_E_yz

static void
cuda2_push_flds_E_yz(struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  int *ldims = patch->ldims;
  for (int iz = -1; iz < ldims[2] + 2; iz++) {
    for (int iy = -1; iy < ldims[1] + 2; iy++) {
      F3(flds, EX, 0,iy,iz) +=
	cny * (F3(flds, HZ, 0,iy,iz) - F3(flds, HZ, 0,iy-1,iz))
	- cnz * (F3(flds, HY, 0,iy,iz) - F3(flds, HY, 0,iy,iz-1))
	- .5 * ppsc->dt * F3(flds, JXI, 0,iy,iz);
    }

    for (int iy = -2; iy < ldims[1] + 2; iy++) {
      F3(flds, EY, 0,iy,iz) +=
	cnz * (F3(flds, HX, 0,iy,iz) - F3(flds, HX, 0,iy,iz-1))
	- .5 * ppsc->dt * F3(flds, JYI, 0,iy,iz);
    }
  }
      
  for (int iz = -2; iz < ldims[2] + 2; iz++) {
    for (int iy = -1; iy < ldims[1] + 2; iy++) {
      F3(flds, EZ, 0,iy,iz) +=
	- cny * (F3(flds, HX, 0,iy,iz) - F3(flds, HX, 0,iy-1,iz))
	- .5 * ppsc->dt * F3(flds, JZI, 0,iy,iz);
    }
  }
}

// ----------------------------------------------------------------------
// cuda2_push_mflds_E_yz

void
cuda2_push_mflds_E_yz(struct psc_mfields *mflds)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    cuda2_push_flds_E_yz(psc_mfields_get_patch(mflds, p));
  }
}

// ----------------------------------------------------------------------
// cuda2_push_flds_H_yz

static void
cuda2_push_flds_H_yz(struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  int *ldims = patch->ldims;
  for (int iz = -1; iz < ldims[2] + 1; iz++) {
    for (int iy = -1; iy < ldims[1] + 1; iy++) {
      F3(flds, HX, 0,iy,iz) -=
	  cny * (F3(flds, EZ, 0,iy+1,iz) - F3(flds, EZ, 0,iy,iz))
	- cnz * (F3(flds, EY, 0,iy,iz+1) - F3(flds, EY, 0,iy,iz));
    }

    for (int iy = -1; iy < ldims[1] + 2; iy++) {
      F3(flds, HY, 0,iy,iz) -=
	  cnz * (F3(flds, EX, 0,iy,iz+1) - F3(flds, EX, 0,iy,iz));
    }
  }
      
  for (int iz = -1; iz < ldims[2] + 2; iz++) {
    for (int iy = -1; iy < ldims[1] + 1; iy++) {
      F3(flds, HZ, 0,iy,iz) -=
	- cny * (F3(flds, EX, 0,iy+1,iz) - F3(flds, EX, 0,iy,iz));
    }
  }
}

// ----------------------------------------------------------------------
// cuda2_push_mflds_H_yz

void
cuda2_push_mflds_H_yz(struct psc_mfields *mflds)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    cuda2_push_flds_H_yz(psc_mfields_get_patch(mflds, p));
  }
}

