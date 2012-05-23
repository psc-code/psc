
#ifndef PSC_GLUE_H
#define PSC_GLUE_H

// Wrappers for Fortran functions

void PIC_push_part_xyz(struct psc *psc, int p, struct psc_particles *pp, struct psc_fields *pf);
void PIC_push_part_xy(struct psc *psc, int p, struct psc_particles *pp, struct psc_fields *pf);
void PIC_push_part_xz(struct psc *psc, int p, struct psc_particles *pp, struct psc_fields *pf);
void PIC_push_part_yz(struct psc *psc, int p, struct psc_particles *pp, struct psc_fields *pf);
void PIC_push_part_z(struct psc *psc, int p, struct psc_particles *pp, struct psc_fields *pf);
void PIC_push_part_z_vay(struct psc *psc, int p, struct psc_particles *pp, struct psc_fields *pf);
void PIC_randomize(struct psc_particles *pp);
void PIC_bin_coll(struct psc_particles *pp);
void PIC_find_cell_indices(struct psc_particles *pp);
void PIC_msa_e(struct psc_fields *pf);
void PIC_msa_h(struct psc_fields *pf);
void PIC_msb_h(struct psc_fields *pf);
void PIC_msb_e(struct psc_fields *pf);
void PIC_pml_msa(struct psc_fields *pf);
void PIC_pml_msb(struct psc_fields *pf);
void PIC_fill_ghosts_h_b(struct psc *psc, int p, struct psc_fields *pf);
void SET_param_pml(struct psc *psc);
void GET_param_domain(void);
void INIT_param_domain(void);
void INIT_param_psc(void);
void INIT_basic(void);
real PSC_p_pulse_z1(real x, real y, real z, real t);
real PSC_s_pulse_z1(real x, real y, real z, real t);
real PSC_p_pulse_z2(real x, real y, real z, real t);
real PSC_s_pulse_z2(real x, real y, real z, real t);

#endif

