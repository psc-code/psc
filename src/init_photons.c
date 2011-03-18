
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
	      
	p->x[0] = xx[0];
	p->x[1] = xx[1];
	p->x[2] = xx[2];
	p->p[0] = photon_np.p[0];
	p->p[1] = photon_np.p[1];
	p->p[2] = photon_np.p[2];
	p->wni  = photon_np.n / photon_np.n_in_cell;
      }
    } psc_foreach_3d_end;

    photons->nr = i;
  }
}

