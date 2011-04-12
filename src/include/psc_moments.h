
#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_moments, struct psc_moments);

void psc_moments_calc_densities(struct psc_moments *moments,
				mfields_base_t *flds, mparticles_base_t *particles,
				mfields_base_t *res);
void psc_moments_calc_v(struct psc_moments *moments,
			mfields_base_t *flds, mparticles_base_t *particles,
			mfields_base_t *res);
void psc_moments_calc_vv(struct psc_moments *moments,
			 mfields_base_t *flds, mparticles_base_t *particles,
			 mfields_base_t *res);

void psc_moments_calc_photon_n(struct psc_moments *moments,
			       mphotons_t *photons, mfields_base_t *res);
