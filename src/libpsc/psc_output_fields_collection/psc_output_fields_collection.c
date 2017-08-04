
#include "psc_output_fields_collection_private.h"

#include "psc_output_fields_private.h"

#include <mrc_io.h>
#include <mrc_params.h>
#include <string.h>
#include <stdlib.h>

// ----------------------------------------------------------------------
// psc_output_fields_collection_set_psc

void
psc_output_fields_collection_set_psc(struct psc_output_fields_collection *output_fields_collection, struct psc *psc)
{
  output_fields_collection->psc = psc;
}

// ----------------------------------------------------------------------
// psc_output_fields_collection_setup

static void
_psc_output_fields_collection_setup(struct psc_output_fields_collection *coll)
{
  char *s = strdup(coll->names), *p;
  char *s_save = s;
  while ((p = strsep(&s, ", "))) {
    struct psc_output_fields *out =
      psc_output_fields_create(psc_output_fields_collection_comm(coll));
    psc_output_fields_set_name(out, p);
    if (strcmp(p, "psc_output_fields") != 0) {
      char s[100];
      sprintf(s, "pfd_%s", p);
      psc_output_fields_set_param_string(out, "pfd", s);
      sprintf(s, "tfd_%s", p);
      psc_output_fields_set_param_string(out, "tfd", s);
    }
    psc_output_fields_set_psc(out, coll->psc);
    psc_output_fields_set_from_options(out);
    psc_output_fields_setup(out);
    psc_output_fields_collection_add_child(coll, (struct mrc_obj *) out);
  }
  free(s);
  free(s_save);
}

// ----------------------------------------------------------------------
// psc_output_fields_collection_write

static void
_psc_output_fields_collection_write(struct psc_output_fields_collection *out, struct mrc_io *io)
{
  mrc_io_write_ref(io, out, "psc", out->psc);
}

// ----------------------------------------------------------------------
// psc_output_fields_collection_add_children_checkpoint

static void
_psc_output_fields_collection_add_children_checkpoint(struct psc_output_fields_collection *coll)
{
  char *s = strdup(coll->names), *p;
  char *s_save = s;
  while ((p = strsep(&s, ", "))) {
    struct psc_output_fields *out =
      psc_output_fields_create(psc_output_fields_collection_comm(coll));
    psc_output_fields_set_name(out, p);
    if (strcmp(p, "psc_output_fields") != 0) {
      char s[100];
      sprintf(s, "pfd_%s", p);
      psc_output_fields_set_param_string(out, "pfd", s);
      sprintf(s, "tfd_%s", p);
      psc_output_fields_set_param_string(out, "tfd", s);
    }
    psc_output_fields_set_psc(out, coll->psc);
    // skip setting-up child; read from io later instead
    psc_output_fields_collection_add_child(coll, (struct mrc_obj *) out);
  }
  free(s);
  free(s_save);
}

// ----------------------------------------------------------------------
// psc_output_fields_collection_read

static void
_psc_output_fields_collection_read(struct psc_output_fields_collection *out, struct mrc_io *io)
{
  out->psc = mrc_io_read_ref(io, out, "psc", psc);
  _psc_output_fields_collection_add_children_checkpoint(out);
  psc_output_fields_collection_read_children(out, io);
}

// ======================================================================
// forward to subclass

void
psc_output_fields_collection_run(struct psc_output_fields_collection *coll,
				 mfields_base_t *flds, struct psc_mparticles *particles)
{
  struct psc_output_fields *out;
  mrc_obj_for_each_child(out, coll, struct psc_output_fields) {
    psc_output_fields_run(out, flds, particles);
  }
}

#define VAR(x) (void *)offsetof(struct psc_output_fields_collection, x)

static struct param psc_output_fields_collection_descr[] = {
  { "names"         , VAR(names)           , PARAM_STRING("psc_output_fields") },
  {},
};
#undef VAR

// ======================================================================
// psc_output_fields_collection class

struct mrc_class_psc_output_fields_collection mrc_class_psc_output_fields_collection = {
  .name             = "psc_output_fields_collection",
  .size             = sizeof(struct psc_output_fields_collection),
  .param_descr      = psc_output_fields_collection_descr,
  .setup            = _psc_output_fields_collection_setup,
  .write            = _psc_output_fields_collection_write,
  .read             = _psc_output_fields_collection_read,
};

