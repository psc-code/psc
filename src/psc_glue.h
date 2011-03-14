
#ifndef PSC_GLUE_H
#define PSC_GLUE_H

// Wrappers for Fortran functions

void PIC_push_part_xyz(int patch, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_xy(int patch, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_xz(int patch, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz(int patch, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_z(int patch, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz_a(int patch, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz_b(int patch, particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_sort(particles_fortran_t *pp);
void PIC_randomize(particles_fortran_t *pp);
void PIC_bin_coll(particles_fortran_t *pp);
void PIC_find_cell_indices(particles_fortran_t *pp);
void PIC_msa(fields_fortran_t *pf);
void PIC_msb(fields_fortran_t *pf);
void PIC_pml_msa(fields_fortran_t *pf);
void PIC_pml_msb(fields_fortran_t *pf);
void OUT_field(void);
void OUT_part(particles_fortran_t *pp);
void CALC_densities(particles_fortran_t *pp, fields_fortran_t *pf);
void SET_param_pml(void);
void PSC_set_patch(int p);
void GET_param_domain(void);
void INIT_param_domain(void);
void INIT_param_psc(void);
f_real **ALLOC_field(void);
void FREE_field(void);
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
void SERV_read_1(int *timestep, int *n_part);
void SERV_read_2(particles_fortran_t *pp, fields_fortran_t *pf);
void SERV_write(particles_fortran_t *pp, fields_fortran_t *pf);

#endif

