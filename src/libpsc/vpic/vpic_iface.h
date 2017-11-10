
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
// vpic_mparticles

struct vpic_mparticles;

struct vpic_mparticles *vpic_mparticles_create();
void vpic_mparticles_ctor_from_simulation(struct vpic_mparticles *vmprts);

// ----------------------------------------------------------------------
// vpic_push_particles

struct vpic_push_particles;

struct vpic_push_particles *vpic_push_particles_create();
void vpic_push_particles_ctor_from_simulation(struct vpic_push_particles *vpushp);
void vpic_push_particles_push_mprts(struct vpic_push_particles *vpushp,
				    struct vpic_mparticles *vmprts,
				    struct vpic_mfields *vmflds);
void vpic_push_particles_prep(struct vpic_push_particles *vpushp,
			      struct vpic_mparticles *mprts,
			      struct vpic_mfields *vmflds);

// ----------------------------------------------------------------------
// vpic_push_fields

void vpic_push_fields_advance_b(struct vpic_mfields *vmflds, double frac);
void vpic_push_fields_advance_e(struct vpic_mfields *vmflds, double frac);

// ----------------------------------------------------------------------
// vpic_marder

struct vpic_marder;

struct vpic_marder *vpic_marder_create();
void vpic_marder_ctor_from_simulation(struct vpic_marder *vmarder);
void vpic_marder_run(struct vpic_marder *vmarder, struct vpic_mfields *vmflds,
		     struct vpic_mparticles *vmprts, int step);

// ----------------------------------------------------------------------
// other (may want an object eventually)

void vpic_sort_run(struct vpic_mparticles *vmprts, int step);
void vpic_collision_run();
void vpic_emitter();
void vpic_current_injection();
void vpic_field_injection();

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

  int status_interval;
};

void vpic_base_init(struct vpic_simulation_info *info);

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
