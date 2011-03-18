
#include <psc.h>
#include <psc_case_private.h>

void
psc_case_init_photons(struct psc_case *_case)
{
  struct psc *psc = _case->psc;

  mphotons_alloc(&psc->mphotons);

  if (!psc_case_ops(_case)->init_photon_np) {
    // if photons aren't initialized, we'll just have zero of them
    return;
  }

  psc_foreach_patch(psc, p) {
    photons_t *photons = &psc->mphotons.p[p];

    int np = 0;
    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      double xx[3] = { CRDX(p, jx), CRDY(p, jy), CRDZ(p, jz) };
      struct psc_photon_np photon_np = {}; // init to all zero
      psc_case_ops(_case)->init_photon_np(_case, xx, &photon_np);
      np += photon_np.n_in_cell;
    } psc_foreach_3d_end;

    photons_alloc(photons, np);

    int i = 0;
    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      double xx[3] = { CRDX(p, jx), CRDY(p, jy), CRDZ(p, jz) };
      struct psc_photon_np photon_np = {}; // init to all zero
      psc_case_ops(_case)->init_photon_np(_case, xx, &photon_np);
	    
      for (int cnt = 0; cnt < photon_np.n_in_cell; cnt++) {
	photon_t *p = photons_get_one(photons, i++);
	      
	p->xi = xx[0];
	p->yi = xx[1];
	p->zi = xx[2];
	p->pxi = photon_np.p[0];
	p->pyi = photon_np.p[1];
	p->pzi = photon_np.p[2];
	p->wni = photon_np.n / photon_np.n_in_cell;
      }
    } psc_foreach_3d_end;

    photons->nr = i;
  }
}

