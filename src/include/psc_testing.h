
#ifndef PSC_TESTING_H
#define PSC_TESTING_H

#include "psc.h"
#include "psc_case.h"

#define assert_equal(x, y, thres) __assert_equal(x, y, #x, #y, thres)

void __assert_equal(double x, double y, const char *xs, const char *ys, double thres);

void psc_testing_dump(struct psc *psc, const char *basename);
void psc_save_particles_ref(struct psc *psc, mparticles_base_t *particles);
void psc_save_fields_ref(struct psc *psc, mfields_base_t *flds);
void psc_check_currents_ref(struct psc *psc, mfields_base_t *flds, double thres);
void psc_check_currents_ref_noghost(struct psc *psc, mfields_base_t *flds, double thres);
void psc_check_fields_ref(struct psc *psc, mfields_base_t *flds, int *m_flds, double thres);
void psc_check_particles_ref(struct psc *psc, mparticles_base_t *particles,
			     double thres, const char *test_str);
void psc_check_particles_sorted(struct psc *psc, mparticles_base_t *particles);

void psc_testing_push_particles(struct psc *psc, const char *s_push_particles);
void psc_testing_save_ref(struct psc *psc);
void psc_testing_push_particles_check(struct psc *psc, double eps_particles, double eps_fields);

struct psc *psc_testing_create_test_yz(const char *s_push_particles, unsigned int mask,
				       char *moments_type);
struct psc_case *psc_create_test_xy(void);
struct psc_case *psc_create_test_xz(void);
struct psc_case *psc_create_test_yz(void);
struct psc_case *psc_create_test_z(void);

// ======================================================================
// psc_test

struct psc_test {
};

void psc_test_create(struct psc *psc);
void psc_test_step(struct psc *psc);

double psc_test_init_field_linear(struct psc *psc, double x[3], int m);
void psc_test_init_npt_rest(struct psc *psc, int kind, double x[3],
			    struct psc_particle_npt *npt);

// ======================================================================
// psc_testing_{init,finalize}

void psc_testing_init(int *argc, char ***argv);
void psc_testing_finalize();

#endif
