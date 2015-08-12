
#include "psc.h"
#include "psc_acc.h"

#include "psc_fields_acc.h"

// ----------------------------------------------------------------------
// acc_push_flds_E_xyz

static void
acc_push_flds_E_xyz(struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_acc_real_t cnx = .5 * ppsc->dt / patch->dx[0];
  fields_acc_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_acc_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  int *ldims = patch->ldims;
  // FIXME, loop limits
  for (int iz = -1; iz < ldims[2] + 2; iz++) {
    for (int iy = -1; iy < ldims[1] + 2; iy++) {
      for (int ix = -1; ix < ldims[0] + 2; ix++) {
      F3_ACC(flds, EX, ix,iy,iz) +=
	+ cny * (F3_ACC(flds, HZ, ix,iy,iz) - F3_ACC(flds, HZ, ix,iy-1,iz))
	- cnz * (F3_ACC(flds, HY, ix,iy,iz) - F3_ACC(flds, HY, ix,iy,iz-1))
	- .5 * ppsc->dt * F3_ACC(flds, JXI, ix,iy,iz);

      F3_ACC(flds, EY, ix,iy,iz) +=
	+ cnz * (F3_ACC(flds, HX, ix,iy,iz) - F3_ACC(flds, HX, ix,iy,iz-1))
	- cnx * (F3_ACC(flds, HZ, ix,iy,iz) - F3_ACC(flds, HZ, ix-1,iy,iz))
	- .5 * ppsc->dt * F3_ACC(flds, JYI, ix,iy,iz);

      F3_ACC(flds, EZ, ix,iy,iz) +=
	+ cnx * (F3_ACC(flds, HY, ix,iy,iz) - F3_ACC(flds, HY, ix-1,iy,iz))
	- cny * (F3_ACC(flds, HX, ix,iy,iz) - F3_ACC(flds, HX, ix,iy-1,iz))
	- .5 * ppsc->dt * F3_ACC(flds, JZI, ix,iy,iz);
      }
    }
  }
}

// ----------------------------------------------------------------------
// acc_push_mflds_E_xyz

void
acc_push_mflds_E_xyz(struct psc_mfields *mflds)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    acc_push_flds_E_xyz(psc_mfields_get_patch(mflds, p));
  }
}

// ----------------------------------------------------------------------
// acc_push_flds_E_yz

static void
acc_push_flds_E_yz(struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_acc_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_acc_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  int *ldims = patch->ldims;
  for (int iz = -1; iz < ldims[2] + 2; iz++) { // FIXME was -1
    for (int iy = -1; iy < ldims[1] + 2; iy++) { // FIXME was -1
      F3_ACC(flds, EX, 0,iy,iz) +=
	cny * (F3_ACC(flds, HZ, 0,iy,iz) - F3_ACC(flds, HZ, 0,iy-1,iz))
	- cnz * (F3_ACC(flds, HY, 0,iy,iz) - F3_ACC(flds, HY, 0,iy,iz-1))
	- .5 * ppsc->dt * F3_ACC(flds, JXI, 0,iy,iz);
    }

    for (int iy = -2; iy < ldims[1] + 2; iy++) {
      F3_ACC(flds, EY, 0,iy,iz) +=
	cnz * (F3_ACC(flds, HX, 0,iy,iz) - F3_ACC(flds, HX, 0,iy,iz-1))
	- .5 * ppsc->dt * F3_ACC(flds, JYI, 0,iy,iz);
    }
  }
      
  for (int iz = -2; iz < ldims[2] + 2; iz++) {
    for (int iy = -1; iy < ldims[1] + 2; iy++) { // FIXME was -1
      F3_ACC(flds, EZ, 0,iy,iz) +=
	- cny * (F3_ACC(flds, HX, 0,iy,iz) - F3_ACC(flds, HX, 0,iy-1,iz))
	- .5 * ppsc->dt * F3_ACC(flds, JZI, 0,iy,iz);
    }
  }
}

// ----------------------------------------------------------------------
// acc_push_mflds_E_yz

void
acc_push_mflds_E_yz(struct psc_mfields *mflds)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    acc_push_flds_E_yz(psc_mfields_get_patch(mflds, p));
  }
}

// ----------------------------------------------------------------------
// acc_push_flds_H_xyz

static void
acc_push_flds_H_xyz(struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_acc_real_t cnx = .5 * ppsc->dt / patch->dx[0];
  fields_acc_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_acc_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  int *ldims = patch->ldims;
  for (int iz = -1; iz < ldims[2] + 1; iz++) {
    for (int iy = -1; iy < ldims[1] + 1; iy++) {
      for (int ix = -1; ix < ldims[0] + 1; ix++) {
	F3_ACC(flds, HX, ix,iy,iz) -=
	  + cny * (F3_ACC(flds, EZ, ix,iy+1,iz) - F3_ACC(flds, EZ, ix,iy,iz))
	  - cnz * (F3_ACC(flds, EY, ix,iy,iz+1) - F3_ACC(flds, EY, ix,iy,iz));
	
	F3_ACC(flds, HY, ix,iy,iz) -=
	  + cnz * (F3_ACC(flds, EX, ix,iy,iz+1) - F3_ACC(flds, EX, ix,iy,iz))
	  - cnx * (F3_ACC(flds, EZ, ix+1,iy,iz) - F3_ACC(flds, EZ, ix,iy,iz));
	
	F3_ACC(flds, HZ, ix,iy,iz) -=
	  + cnx * (F3_ACC(flds, EY, ix+1,iy,iz) - F3_ACC(flds, EY, ix,iy,iz))
	  - cny * (F3_ACC(flds, EX, ix,iy+1,iz) - F3_ACC(flds, EX, ix,iy,iz));
      }
    }
  }
}

// ----------------------------------------------------------------------
// acc_push_mflds_H_xyz

void
acc_push_mflds_H_xyz(struct psc_mfields *mflds)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    acc_push_flds_H_xyz(psc_mfields_get_patch(mflds, p));
  }
}

// ----------------------------------------------------------------------
// acc_push_flds_H_yz

static void
acc_push_flds_H_yz(struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_acc_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_acc_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  int *ldims = patch->ldims;
  for (int iz = -1; iz < ldims[2] + 1; iz++) { // FIXME was -/+ 1
    for (int iy = -1; iy < ldims[1] + 1; iy++) { // FIXME was -/+ 1
      F3_ACC(flds, HX, 0,iy,iz) -=
	  cny * (F3_ACC(flds, EZ, 0,iy+1,iz) - F3_ACC(flds, EZ, 0,iy,iz))
	- cnz * (F3_ACC(flds, EY, 0,iy,iz+1) - F3_ACC(flds, EY, 0,iy,iz));
    }

    for (int iy = -1; iy < ldims[1] + 2; iy++) { // FIXME was -1 / +2
      F3_ACC(flds, HY, 0,iy,iz) -=
	  cnz * (F3_ACC(flds, EX, 0,iy,iz+1) - F3_ACC(flds, EX, 0,iy,iz));
    }
  }
      
  for (int iz = -1; iz < ldims[2] + 2; iz++) { // FIXME was -1 / + 2
    for (int iy = -1; iy < ldims[1] + 1; iy++) { // FIXME was -/+ 1
      F3_ACC(flds, HZ, 0,iy,iz) -=
	- cny * (F3_ACC(flds, EX, 0,iy+1,iz) - F3_ACC(flds, EX, 0,iy,iz));
    }
  }
}

// ----------------------------------------------------------------------
// acc_push_mflds_H_yz

void
acc_push_mflds_H_yz(struct psc_mfields *mflds)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    acc_push_flds_H_yz(psc_mfields_get_patch(mflds, p));
  }
}

