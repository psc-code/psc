
#ifndef PSC_TESTING_H
#define PSC_TESTING_H

#include "psc.h"

#define assert_equal(x, y, thres) __assert_equal(x, y, #x, #y, thres)

void __assert_equal(double x, double y, const char *xs, const char *ys, double thres);

void psc_save_particles_ref();
void psc_save_fields_ref();
void psc_check_currents_ref(double thres);
void psc_check_currents_ref_noghost(double thres);
void psc_check_fields_ref(int *flds, double thres);
void psc_check_particles_ref(double thres);
void psc_check_particles_sorted();

void psc_create_test_xz(struct psc_mod_config *conf);
void psc_create_test_yz(struct psc_mod_config *conf);

#endif
