
#ifndef VPIC_IFACE_H
#define VPIC_IFACE_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0 // hack to fix indentation
}
#endif

// ----------------------------------------------------------------------
// vpic_mfields

struct vpic_mfields;

struct vpic_mfields *vpic_mfields_create();
void vpic_mfields_ctor_from_simulation(struct vpic_mfields *vmflds);

// ----------------------------------------------------------------------

struct vpic_info {
  double dx, dy, dz;
  double dt;
  double c;
  double eps0;
};

struct vpic_simulation_info {
  int num_step;
  int clean_div_e_interval;
  int clean_div_b_interval;
  int sync_shared_interval;
};

void vpic_base_init(struct vpic_simulation_info *info);
void vpic_base_integrate();

bool vpic_done();
void vpic_performance_sort();
void vpic_clear_accumulator_array();
void vpic_collisions();
void vpic_advance_p();
void vpic_emitter();
void vpic_reduce_accumulator_array();
void vpic_boundary_p(struct vpic_mfields *vmflds);
void vpic_calc_jf(struct vpic_mfields *vmflds);
void vpic_current_injection();
void vpic_advance_b(struct vpic_mfields *vmflds, double frac);
void vpic_advance_e(struct vpic_mfields *vmflds, double frac);
void vpic_field_injection();
void vpic_clean_div_e();
void vpic_clean_div_b();
void vpic_sync_faces();
void vpic_load_interpolator_array();
void vpic_print_status();
void vpic_diagnostics();
void vpic_inc_step(int step);


#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
