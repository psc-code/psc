
#pragma once

// ======================================================================
// SetupFields

template<typename MF>
struct SetupFields
{
  using Mfields = MF;
  using Fields = typename Mfields::fields_t;
  using real_t = typename Mfields::real_t;
  
  static void set_ic(struct psc *psc)
  {
    double (*init_field)(struct psc *psc, double x[3], int m);
    init_field = psc_ops(psc)->init_field;
    if (!init_field)
      return;

    auto mflds_base = PscMfieldsBase{psc->flds};
    mfields_t mf = mflds_base.get_as<mfields_t>(0, 0);

    // FIXME, do we need the ghost points?
    psc_foreach_patch(psc, p) {
      Fields F(mf[p]);

      psc_foreach_3d_g(psc, p, jx, jy, jz) {
	double dx = psc->grid().dx[0], dy = psc->grid().dx[1], dz = psc->grid().dx[2];
	double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);

	double ncc[3] = { xx        , yy + .5*dy, zz + .5*dz };
	double cnc[3] = { xx + .5*dx, yy        , zz + .5*dz };
	double ccn[3] = { xx + .5*dx, yy + .5*dy, zz         };

	double cnn[3] = { xx + .5*dx, yy        , zz         };
	double ncn[3] = { xx        , yy + .5*dy, zz         };
	double nnc[3] = { xx        , yy        , zz + .5*dz };
 
	F(HX, jx,jy,jz) += init_field(psc, ncc, HX);
	F(HY, jx,jy,jz) += init_field(psc, cnc, HY);
	F(HZ, jx,jy,jz) += init_field(psc, ccn, HZ);

	F(EX, jx,jy,jz) += init_field(psc, cnn, EX);
	F(EY, jx,jy,jz) += init_field(psc, ncn, EY);
	F(EZ, jx,jy,jz) += init_field(psc, nnc, EZ);

	F(JXI, jx,jy,jz) += init_field(psc, cnn, JXI);
	F(JYI, jx,jy,jz) += init_field(psc, ncn, JYI);
	F(JZI, jx,jy,jz) += init_field(psc, nnc, JZI);

      } foreach_3d_g_end;
    }
    mf.put_as(mflds_base, JXI, HX + 3);
  }
};
