
#ifndef PSC_TESTING_H
#define PSC_TESTING_H

#include "psc.h"

#define assert_equal(x, y, thres) __assert_equal(x, y, #x, #y, thres)

void __assert_equal(double x, double y, const char *xs, const char *ys, double thres);

void psc_save_particles_ref(mparticles_base_t *particles);
void psc_save_fields_ref(mfields_base_t *flds);
void psc_check_currents_ref(mfields_base_t *flds, double thres);
void psc_check_currents_ref_noghost(mfields_base_t *flds, double thres);
void psc_check_fields_ref(mfields_base_t *flds, int *m_flds, double thres);
void psc_check_particles_ref(mparticles_base_t *particles, double thres, const char *test_str);
void psc_check_particles_sorted(mparticles_base_t *particles);

void psc_create_test_xy(struct psc_mod_config *conf);
void psc_create_test_xz(struct psc_mod_config *conf);
void psc_create_test_yz(struct psc_mod_config *conf);

#endif
