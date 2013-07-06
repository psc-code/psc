
#include "ggcm_mhd_crds_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds_gen.h"

#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>

static const char *crdname[NR_CRDS] = {
  [FX1] = "FX1",
  [FD1] = "FD1",
  [FX2] = "FX2",
  [BD1] = "BD1",
  [BD2] = "BD2",
  [BD3] = "BD3",
  [BD4] = "BD4",
};

// ======================================================================
// ggcm_mhd_crds class

#define ggcm_mhd_crds_ops(crds) ((struct ggcm_mhd_crds_ops *)((crds)->obj.ops))

// ----------------------------------------------------------------------
// ggcm_mhd_crds_create

static void
_ggcm_mhd_crds_create(struct ggcm_mhd_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    crds->f1[d] = mrc_f1_create(MPI_COMM_SELF);
    char s[10]; sprintf(s, "f1[%d]", d);
    mrc_f1_set_name(crds->f1[d], s);
    mrc_f1_set_nr_comps(crds->f1[d], NR_CRDS);
    mrc_f1_set_param_int(crds->f1[d], "dim", d);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_setup

static void
_ggcm_mhd_crds_setup(struct ggcm_mhd_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    mrc_f1_setup(crds->f1[d]);
    for (int m = 0; m < NR_CRDS; m++) {
      mrc_f1_set_comp_name(crds->f1[d], m, crdname[m]);
    }
  }

  // set up values for Fortran coordinate arrays
  assert(crds->crds_gen);
  ggcm_mhd_crds_gen_run(crds->crds_gen, crds);

  // set up mrc_crds
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(crds->domain, 0, &info);

  struct mrc_crds *mrc_crds = mrc_domain_get_crds(crds->domain);
  mrc_crds_set_values(mrc_crds,
		      ggcm_mhd_crds_get_crd(crds, 0, FX1), info.ldims[0],
		      ggcm_mhd_crds_get_crd(crds, 1, FX1), info.ldims[1],
		      ggcm_mhd_crds_get_crd(crds, 2, FX1), info.ldims[2]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_destroy

static void
_ggcm_mhd_crds_destroy(struct ggcm_mhd_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    mrc_f1_destroy(crds->f1[d]);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_read

static void
_ggcm_mhd_crds_read(struct ggcm_mhd_crds *crds, struct mrc_io *io)
{
  ggcm_mhd_crds_read_member_objs(crds, io);
  crds->f1[0] = mrc_io_read_ref(io, crds, "f1[0]", mrc_f1);
  crds->f1[1] = mrc_io_read_ref(io, crds, "f1[1]", mrc_f1);
  crds->f1[2] = mrc_io_read_ref(io, crds, "f1[2]", mrc_f1);
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_write

static void
_ggcm_mhd_crds_write(struct ggcm_mhd_crds *crds, struct mrc_io *io)
{
  mrc_io_write_ref(io, crds, "f1[0]", crds->f1[0]);
  mrc_io_write_ref(io, crds, "f1[1]", crds->f1[1]);
  mrc_io_write_ref(io, crds, "f1[2]", crds->f1[2]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_get_crd

float *
ggcm_mhd_crds_get_crd(struct ggcm_mhd_crds *crds, int d, int m)
{
  return &MRC_F1(crds->f1[d], m, 0);
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_init

static void
ggcm_mhd_crds_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_crds, &ggcm_mhd_crds_c_ops);
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds class description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_crds, x)
static struct param ggcm_mhd_crds_descr[] = {
  { "domain"          , VAR(domain)         , PARAM_OBJ(mrc_domain)          },

  { "crds_gen"        , VAR(crds_gen)       , MRC_VAR_OBJ(ggcm_mhd_crds_gen) },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_crds class

struct mrc_class_ggcm_mhd_crds mrc_class_ggcm_mhd_crds = {
  .name             = "ggcm_mhd_crds",
  .size             = sizeof(struct ggcm_mhd_crds),
  .param_descr      = ggcm_mhd_crds_descr,
  .init             = ggcm_mhd_crds_init,
  .create           = _ggcm_mhd_crds_create,
  .setup            = _ggcm_mhd_crds_setup,
  .destroy          = _ggcm_mhd_crds_destroy,
  .read             = _ggcm_mhd_crds_read,
  .write            = _ggcm_mhd_crds_write,
};

