
#pragma once

// ======================================================================
// SetupFields

template<typename MF>
struct SetupFields
{
  using Mfields = MF;
  using Fields = typename Mfields::fields_t;
  using real_t = typename Mfields::real_t;

  template<typename FUNC>
  static void set(Mfields& mf, FUNC func)
  {
    psc_foreach_patch(ppsc, p) {
      Fields F(mf[p]);

      // FIXME, do we need the ghost points?
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	struct psc* psc = ppsc;
	double dx = mf.grid().dx[0], dy = mf.grid().dx[1], dz = mf.grid().dx[2];
	double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);

	double ncc[3] = { xx        , yy + .5*dy, zz + .5*dz };
	double cnc[3] = { xx + .5*dx, yy        , zz + .5*dz };
	double ccn[3] = { xx + .5*dx, yy + .5*dy, zz         };

	double cnn[3] = { xx + .5*dx, yy        , zz         };
	double ncn[3] = { xx        , yy + .5*dy, zz         };
	double nnc[3] = { xx        , yy        , zz + .5*dz };

	F(HX, jx,jy,jz) += func(HX, ncc);
	F(HY, jx,jy,jz) += func(HY, cnc);
	F(HZ, jx,jy,jz) += func(HZ, ccn);

	F(EX, jx,jy,jz) += func(EX, cnn);
	F(EY, jx,jy,jz) += func(EY, ncn);
	F(EZ, jx,jy,jz) += func(EZ, nnc);

	F(JXI, jx,jy,jz) += func(JXI, cnn);
	F(JYI, jx,jy,jz) += func(JYI, ncn);
	F(JZI, jx,jy,jz) += func(JZI, nnc);
      } foreach_3d_g_end;
    }
  }
  
  static void set_ic(struct psc *psc)
  {
    double (*init_field)(struct psc *psc, double x[3], int m);
    init_field = psc_ops(psc)->init_field;
    if (!init_field)
      return;

    auto mflds_base = PscMfieldsBase{psc->flds};
    mfields_t mflds = mflds_base.get_as<mfields_t>(0, 0);

    set(*mflds.sub(), [&](int m, real_t xx[3]) {
	return init_field(psc, xx, m);
      });

    mflds.put_as(mflds_base, JXI, HX + 3);
  }
};
