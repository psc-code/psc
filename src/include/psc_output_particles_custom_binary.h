/*
 *  psc_output_particles_c.c
 *  psc
 *
 *  Created by Nils Moschuering on 26.09.11.
 *  Copyright 2011 LMU. All rights reserved.
 *
 */

#ifndef PSC_OUTPUT_PARTICLES_CUSTOM_BINARY_H
#define PSC_OUTPUT_PARTICLES_CUSTOM_BINARY_H

#include "psc_output_particles_private.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_c.h"

#include <stdbool.h>

struct psc_output_particles_custom_binary {
  char *data_dir;
  struct mrc_io *io;
  bool (*filter_func)(particle_t *part);
  int first;
  int step;
  int next;
  bool write_q;
  bool write_m;
  bool write_x;
  bool write_y;
  bool write_z;
  bool write_px;
  bool write_py;
  bool write_pz;
};

void
psc_output_particles_custom_binary_setfilter(struct psc_output_particles *out, bool (*new_filter_func)(particle_t *part));

static void write_particles_to_file(struct psc_output_particles_custom_binary *out_c, mparticles_base_t *particles);

static void create_output_file(struct psc_output_particles_custom_binary *out_c);

#endif
