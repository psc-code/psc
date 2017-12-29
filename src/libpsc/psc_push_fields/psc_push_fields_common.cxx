
#include "psc.h"

#include "psc_push_fields_impl.hxx"

struct params_push_fields {
  fields_real_t dth;
  fields_real_t cnx;
  fields_real_t cny;
  fields_real_t cnz;
  int ldims[3];
};

static struct params_push_fields prm;

// ----------------------------------------------------------------------
// params_push_fields_set

static void
params_push_fields_set(struct psc *psc, double dt_fac)
{
  for (int p = 1; p < psc->nr_patches; p++) {
    for (int d = 0; d < 3; d++) {
      assert(psc->patch[0].ldims[d] == psc->patch[p].ldims[d]);
      assert(psc->patch[0].dx[d] == psc->patch[p].dx[d]);
    }
  }

  prm.dth = dt_fac * psc->dt;

  prm.cnx = prm.dth / psc->patch[0].dx[0];
  prm.cny = prm.dth / psc->patch[0].dx[1];
  prm.cnz = prm.dth / psc->patch[0].dx[2];

  if (psc->domain.gdims[0] == 1) {
    prm.cnx = 0.;
  }
  if (psc->domain.gdims[1] == 1) {
    prm.cny = 0.;
  }
  if (psc->domain.gdims[2] == 1) {
    prm.cnz = 0.;
  }

  for (int d = 0; d < 3; d++) {
    prm.ldims[d] = psc->patch[0].ldims[d];
  }
}

#define MAKE_PUSH_E(func)						\
									\
static void								\
func(struct psc_push_fields *push, fields_t flds)			\
{									\
  foreach_3d(ppsc, 0, i,j,k, 1, 2) {					\
    F3(flds, EX, i,j,k) +=						\
      prm.cny * (F3(flds, HZ, i,j,k) - F3(flds, HZ, i,j-1,k)) -		\
      prm.cnz * (F3(flds, HY, i,j,k) - F3(flds, HY, i,j,k-1)) -		\
      prm.dth * F3(flds, JXI, i,j,k);					\
    									\
    F3(flds, EY, i,j,k) +=						\
      prm.cnz * (F3(flds, HX, i,j,k) - F3(flds, HX, i,j,k-1)) -		\
      prm.cnx * (F3(flds, HZ, i,j,k) - F3(flds, HZ, i-1,j,k)) -		\
      prm.dth * F3(flds, JYI, i,j,k);					\
    									\
    F3(flds, EZ, i,j,k) +=						\
      prm.cnx * (F3(flds, HY, i,j,k) - F3(flds, HY, i-1,j,k)) -		\
      prm.cny * (F3(flds, HX, i,j,k) - F3(flds, HX, i,j-1,k)) -		\
      prm.dth * F3(flds, JZI, i,j,k);					\
  } foreach_3d_end;							\
}

static void
psc_push_fields_sub_push_E_yz(struct psc_push_fields *push, fields_t flds)
{
  for (int iz = -1; iz < prm.ldims[2] + 2; iz++) {
    for (int iy = -1; iy < prm.ldims[1] + 2; iy++) {
      _F3(flds, EX, 0,iy,iz) +=
	prm.cny * (_F3(flds, HZ, 0,iy,iz) - _F3(flds, HZ, 0,iy-1,iz)) -
	prm.cnz * (_F3(flds, HY, 0,iy,iz) - _F3(flds, HY, 0,iy,iz-1)) -
	prm.dth * _F3(flds, JXI, 0,iy,iz);
    }

    for (int iy = -2; iy < prm.ldims[1] + 2; iy++) {
      _F3(flds, EY, 0,iy,iz) +=
	prm.cnz * (_F3(flds, HX, 0,iy,iz) - _F3(flds, HX, 0,iy,iz-1)) -
	0 -
	prm.dth * _F3(flds, JYI, 0,iy,iz);
    }
  }
      
  for (int iz = -2; iz < prm.ldims[2] + 2; iz++) {
    for (int iy = -1; iy < prm.ldims[1] + 2; iy++) {
      _F3(flds, EZ, 0,iy,iz) +=
	0 -
	prm.cny * (_F3(flds, HX, 0,iy,iz) - _F3(flds, HX, 0,iy-1,iz)) -
	prm.dth * _F3(flds, JZI, 0,iy,iz);
    }
  }
}

#define MAKE_PUSH_H(func)						\
									\
static void								\
func(struct psc_push_fields *push, fields_t flds)			\
{									\
  foreach_3d(ppsc, 0, i,j,k, 2, 1) {					\
    F3(flds, HX, i,j,k) -=						\
      prm.cny * (F3(flds, EZ, i,j+1,k) - F3(flds, EZ, i,j,k)) -		\
      prm.cnz * (F3(flds, EY, i,j,k+1) - F3(flds, EY, i,j,k));		\
									\
    F3(flds, HY, i,j,k) -=						\
      prm.cnz * (F3(flds, EX, i,j,k+1) - F3(flds, EX, i,j,k)) -		\
      prm.cnx * (F3(flds, EZ, i+1,j,k) - F3(flds, EZ, i,j,k));		\
									\
    F3(flds, HZ, i,j,k) -=						\
      prm.cnx * (F3(flds, EY, i+1,j,k) - F3(flds, EY, i,j,k)) -		\
      prm.cny * (F3(flds, EX, i,j+1,k) - F3(flds, EX, i,j,k));		\
  } foreach_3d_end;							\
}

#define F3(flds, m, i,j,k) _F3(flds, m, i,j,k)
MAKE_PUSH_E(psc_push_fields_sub_push_E_xyz)
MAKE_PUSH_H(psc_push_fields_sub_push_H_xyz)
#undef F3

static void
psc_push_fields_sub_push_H_yz(struct psc_push_fields *push, fields_t flds)
{
  for (int iz = -1; iz < prm.ldims[2] + 1; iz++) {
    for (int iy = -1; iy < prm.ldims[1] + 1; iy++) {
      _F3(flds, HX, 0,iy,iz) -=
	prm.cny * (_F3(flds, EZ, 0,iy+1,iz) - _F3(flds, EZ, 0,iy,iz)) -
	prm.cnz * (_F3(flds, EY, 0,iy,iz+1) - _F3(flds, EY, 0,iy,iz));
    }

    for (int iy = -1; iy < prm.ldims[1] + 2; iy++) {
      _F3(flds, HY, 0,iy,iz) -=
	prm.cnz * (_F3(flds, EX, 0,iy,iz+1) - _F3(flds, EX, 0,iy,iz)) -
	0;
    }
  }
      
  for (int iz = -1; iz < prm.ldims[2] + 2; iz++) {
    for (int iy = -1; iy < prm.ldims[1] + 1; iy++) {
      _F3(flds, HZ, 0,iy,iz) -=
	0 -
	prm.cny * (_F3(flds, EX, 0,iy+1,iz) - _F3(flds, EX, 0,iy,iz));
    }
  }
}

// E-field propagation E^(n)    , H^(n), j^(n) 
//                  -> E^(n+0.5), H^(n), j^(n)
// Ex^{n}[-.5:+.5][-1:1][-1:1] -> Ex^{n+.5}[-.5:+.5][-1:1][-1:1]
// using Hx^{n}[-1:1][-1.5:1.5][-1.5:1.5]
//       jx^{n+1}[-.5:.5][-1:1][-1:1]

static void
psc_push_fields_sub_push_mflds_E(struct psc_push_fields *push, struct psc_mfields *mflds_base,
				 double dt_fac)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, JXI, HX + 3);

  params_push_fields_set(ppsc, dt_fac);
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    int *gdims = ppsc->domain.gdims;
    if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_sub_push_E_xyz(push, flds);
    } else if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_sub_push_E_yz(push, flds);
    } else if (gdims[0] > 1 && gdims[1] == 1 && gdims[2] > 1) {
#ifdef PSC_FIELDS_AS_SINGLE
      using Fields = Fields3d<fields_single_real_t, fields_single_t, DIM_XZ>;
      psc_push_fields_push_E<Fields>(push, flds, ppsc, dt_fac);
#elif PSC_FIELDS_AS_C
      using Fields = Fields3d<fields_c_real_t, fields_c_t, DIM_XZ>;
      psc_push_fields_push_E<Fields>(push, flds, ppsc, dt_fac);
#else
      psc_push_fields_sub_push_E_xz(push, flds);
#endif
    } else {
      assert(0);
    }
  }

  psc_mfields_put_as(mflds, mflds_base, EX, EX + 3);
}

// B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
//                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)
// Hx^{n}[:][-.5:.5][-.5:.5] -> Hx^{n+.5}[:][-.5:.5][-.5:.5]
// using Ex^{n+.5}[-.5:+.5][-1:1][-1:1]

static void
psc_push_fields_sub_push_mflds_H(struct psc_push_fields *push, struct psc_mfields *mflds_base,
				 double dt_fac)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, HX + 3);

  params_push_fields_set(ppsc, dt_fac);
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    int *gdims = ppsc->domain.gdims;
    if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_sub_push_H_xyz(push, flds);
    } else if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_sub_push_H_yz(push, flds);
    } else if (gdims[0] > 1 && gdims[1] == 1 && gdims[2] > 1) {
#ifdef PSC_FIELDS_AS_SINGLE
      using Fields = Fields3d<fields_single_real_t, fields_single_t, DIM_XZ>;
      psc_push_fields_push_H<Fields>(push, flds, ppsc, dt_fac);
#elif PSC_FIELDS_AS_C
      using Fields = Fields3d<fields_c_real_t, fields_c_t, DIM_XZ>;
      psc_push_fields_push_H<Fields>(push, flds, ppsc, dt_fac);
#else
      psc_push_fields_sub_push_H_xz(push, flds);
#endif
    } else {
      assert(0);
    }
  }

  psc_mfields_put_as(mflds, mflds_base, HX, HX + 3);
}

