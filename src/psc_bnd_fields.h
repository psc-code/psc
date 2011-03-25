
#ifndef PSC_BND_FIELDS_H
#define PSC_BND_FIELDS_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_bnd_fields, struct psc_bnd_fields);

struct psc_pulse *psc_bnd_fields_get_pulse_x1(struct psc_bnd_fields *bnd);
struct psc_pulse *psc_bnd_fields_get_pulse_x2(struct psc_bnd_fields *bnd);
struct psc_pulse *psc_bnd_fields_get_pulse_y1(struct psc_bnd_fields *bnd);
struct psc_pulse *psc_bnd_fields_get_pulse_y2(struct psc_bnd_fields *bnd);
struct psc_pulse *psc_bnd_fields_get_pulse_z1(struct psc_bnd_fields *bnd);
struct psc_pulse *psc_bnd_fields_get_pulse_z2(struct psc_bnd_fields *bnd);
void psc_bnd_fields_setup_fields(struct psc_bnd_fields *bnd, mfields_base_t *mflds);
void psc_bnd_fields_fill_ghosts_b_H(struct psc_bnd_fields *bnd, mfields_base_t *mflds);

#endif
