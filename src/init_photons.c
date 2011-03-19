
#include <psc.h>
#include <psc_case_private.h>

#include <math.h>

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
	      
	float ran1, ran2, ran3, ran4, ran5, ran6;
	do {
	  ran1 = random() / ((float) RAND_MAX + 1);
	  ran2 = random() / ((float) RAND_MAX + 1);
	  ran3 = random() / ((float) RAND_MAX + 1);
	  ran4 = random() / ((float) RAND_MAX + 1);
	  ran5 = random() / ((float) RAND_MAX + 1);
	  ran6 = random() / ((float) RAND_MAX + 1);
	} while (ran1 >= 1.f || ran2 >= 1.f || ran3 >= 1.f ||
		 ran4 >= 1.f || ran5 >= 1.f || ran6 >= 1.f);
	
	float kx =
	  sqrtf(-2.f*photon_np.sigma_k[0]*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2)
	  + photon_np.k[0];
	float ky =
	  sqrtf(-2.f*photon_np.sigma_k[1]*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4)
	  + photon_np.k[1];
	float kz =
	  sqrtf(-2.f*photon_np.sigma_k[2]*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6)
	  + photon_np.k[2];

	p->x[0] = xx[0];
	p->x[1] = xx[1];
	p->x[2] = xx[2];
	p->p[0] = kx;
	p->p[1] = ky;
	p->p[2] = kz;
	p->wni  = photon_np.n / photon_np.n_in_cell;
      }
    } psc_foreach_3d_end;

    photons->nr = i;
  }
}

