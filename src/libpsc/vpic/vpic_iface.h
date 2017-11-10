
#ifndef VPIC_IFACE_H
#define VPIC_IFACE_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0 // hack to fix indentation
}
#endif

struct vpic_info {
  double dx, dy, dz;
  double dt;
  double c;
  double eps0;
};

void vpic_base_init(struct vpic_info *info);
void vpic_base_integrate();

bool vpic_done();
void vpic_performance_sort();
void vpic_clear_accumulator_array();
void vpic_collisions();
void vpic_advance_p();
void vpic_emitter();
void vpic_reduce_accumulator_array();
void vpic_boundary_p();
void vpic_calc_jf();
void vpic_current_injection();
void vpic_advance_b(double frac);
void vpic_advance_e(double frac);
void vpic_field_injection();
void vpic_clean_div_e();
void vpic_clean_div_b();
void vpic_sync_faces();
void vpic_load_interpolator_array();
void vpic_print_status();
void vpic_diagnostics();
void vpic_step();


#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
