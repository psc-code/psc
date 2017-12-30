
#ifndef PSC_TESTING_H
#define PSC_TESTING_H

#include "psc.h"

BEGIN_C_DECLS

#define assert_equal(x, y, thres) __assert_equal(x, y, #x, #y, thres)

void __assert_equal(double x, double y, const char *xs, const char *ys, double thres);

void psc_testing_dump(struct psc *psc, const char *basename);
void psc_save_particles_ref(struct psc *psc, struct psc_mparticles *particles);
void psc_save_fields_ref(struct psc *psc, struct psc_mfields *flds);
void psc_check_currents_ref(struct psc *psc, struct psc_mfields *flds, double thres, int sw);
void psc_check_fields_ref(struct psc *psc, struct psc_mfields *flds, int *m_flds, double thres);
void psc_check_particles_ref(struct psc *psc, struct psc_mparticles *particles,
			     double thres, const char *test_str);
void psc_check_particles_sorted(struct psc *psc, struct psc_mparticles *particles);
void psc_testing_check_densities_ref(struct psc *psc, struct psc_mparticles *particles,
				     double eps);

void psc_testing_push_particles(struct psc *psc, const char *s_push_particles);
void psc_testing_save_ref(struct psc *psc);
void psc_testing_push_particles_check(struct psc *psc, double eps_particles, double eps_fields);

// test with linear E, B, particles at rest
extern struct psc_ops psc_test_ops_1;

struct psc *psc_testing_create_test_yz(const char *s_push_particles, unsigned int mask);
struct psc *psc_testing_create_test_xz();

// ======================================================================
// psc_test

struct psc_test {
  int dummy;
};

void psc_test_create(struct psc *psc);
void psc_test_step(struct psc *psc);

double psc_test_init_field_linear(struct psc *psc, double x[3], int m);
void psc_test_init_npt_rest(struct psc *psc, int kind, double x[3],
			    struct psc_particle_npt *npt);

// ======================================================================
// psc_testing_{init,finalize}

void psc_testing_init(int *argc, char ***argv);
void psc_testing_finalize(void);

extern bool opt_testing_dump;

END_C_DECLS

#endif
