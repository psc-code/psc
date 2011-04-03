
#ifndef PSC_TESTING_H
#define PSC_TESTING_H

#include "psc.h"
#include "psc_case.h"

#define assert_equal(x, y, thres) __assert_equal(x, y, #x, #y, thres)

void __assert_equal(double x, double y, const char *xs, const char *ys, double thres);

void psc_save_particles_ref(struct psc *psc, mparticles_base_t *particles);
void psc_save_fields_ref(struct psc *psc, mfields_base_t *flds);
void psc_check_currents_ref(struct psc *psc, mfields_base_t *flds, double thres);
void psc_check_currents_ref_noghost(struct psc *psc, mfields_base_t *flds, double thres);
void psc_check_fields_ref(struct psc *psc, mfields_base_t *flds, int *m_flds, double thres);
void psc_check_particles_ref(struct psc *psc, mparticles_base_t *particles,
			     double thres, const char *test_str);
void psc_check_particles_sorted(struct psc *psc, mparticles_base_t *particles);

struct psc_case *psc_create_test_xz(void);
struct psc_case *psc_create_test_yz(void);

#endif
