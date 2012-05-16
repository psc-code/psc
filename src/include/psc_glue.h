
#ifndef PSC_GLUE_H
#define PSC_GLUE_H

// Wrappers for Fortran functions

void PIC_push_part_xyz(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_xy(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_xz(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_z(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz_a(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz_b(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_z_vay(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_sort(particles_fortran_t *pp);
void PIC_randomize(particles_fortran_t *pp);
void PIC_bin_coll(particles_fortran_t *pp);
void PIC_find_cell_indices(particles_fortran_t *pp);
void PIC_msa_e(fields_fortran_t *pf);
void PIC_msa_h(fields_fortran_t *pf);
void PIC_msb_h(fields_fortran_t *pf);
void PIC_msb_e(fields_fortran_t *pf);
void PIC_pml_msa(fields_fortran_t *pf);
void PIC_pml_msb(fields_fortran_t *pf);
void PIC_fill_ghosts_h_b(struct psc *psc, int p, fields_fortran_t *pf);
void SET_param_pml(struct psc *psc);
void GET_param_domain(void);
void INIT_param_domain(void);
void INIT_param_psc(void);
void INIT_basic(void);
real PSC_p_pulse_z1(real x, real y, real z, real t);
real PSC_s_pulse_z1(real x, real y, real z, real t);
real PSC_p_pulse_z2(real x, real y, real z, real t);
real PSC_s_pulse_z2(real x, real y, real z, real t);

void PIC_fax(fields_fortran_t *pf, int m);
void PIC_fay(fields_fortran_t *pf, int m);
void PIC_faz(fields_fortran_t *pf, int m);
void PIC_fex(fields_fortran_t *pf, int m);
void PIC_fey(fields_fortran_t *pf, int m);
void PIC_fez(fields_fortran_t *pf, int m);
void PIC_pex(particles_fortran_t *pp);
void PIC_pey(particles_fortran_t *pp);
void PIC_pez(particles_fortran_t *pp);

#endif

