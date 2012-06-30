
#include <mrc_a3.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_params.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>


// ======================================================================
// mrc_a3

static void
_mrc_a3_destroy(struct mrc_a3 *a3)
{
  for (int p = 0; p < a3->nr_patches; p++) {
    struct mrc_a3_patch *a3p = &a3->patches[p];
    free(a3p->arr);
  }
  free(a3->patches);

  if (a3->name) {
    for (int m = 0; m < a3->nr_comp; m++) {
      free(a3->name[m]);
    }
    free(a3->name);
  }
}

static void
_mrc_a3_setup(struct mrc_a3 *a3)
{
  a3->name = calloc(a3->nr_comp, sizeof(*a3->name));

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(a3->domain, &nr_patches);

  a3->nr_patches = nr_patches;
  a3->patches = calloc(nr_patches, sizeof(*a3->patches));
  for (int p = 0; p < nr_patches; p++) {
    struct mrc_a3_patch *a3p = &a3->patches[p];
    for (int d = 0; d < 3; d++) {
      a3p->ib[d] = -a3->sw;
      a3p->im[d] = patches[p].ldims[d] + 2 * a3->sw;
    }
    int len = a3p->im[0] * a3p->im[1] * a3p->im[2] * a3->nr_comp;
    a3p->arr = calloc(len, sizeof(*a3p->arr));
  }
}

static void
_mrc_a3_view(struct mrc_a3 *a3)
{
#if 1
  int rank, size;
  MPI_Comm_rank(mrc_a3_comm(a3), &rank);
  MPI_Comm_size(mrc_a3_comm(a3), &size);

  for (int r = 0; r < size; r++) {
    if (r == rank) {
      mrc_a3_foreach_patch(a3, p) {
	struct mrc_a3_patch *a3p = mrc_a3_patch_get(a3, p);
	mprintf("patch %d: ib = %dx%dx%d im = %dx%dx%d\n", p,
		a3p->ib[0], a3p->ib[1], a3p->ib[2],
		a3p->im[0], a3p->im[1], a3p->im[2]);
      }
    }
    MPI_Barrier(mrc_a3_comm(a3));
  }
#endif
}

void
mrc_a3_set_comp_name(struct mrc_a3 *a3, int m, const char *name)
{
  assert(m < a3->nr_comp);
  free(a3->name[m]);
  a3->name[m] = name ? strdup(name) : NULL;
}

const char *
mrc_a3_comp_name(struct mrc_a3 *a3, int m)
{
  assert(m < a3->nr_comp);
  return a3->name[m];
}

static void
_mrc_a3_write(struct mrc_a3 *a3, struct mrc_io *io)
{
  mrc_io_write_obj_ref(io, mrc_a3_name(a3), "domain",
		       (struct mrc_obj *) a3->domain);
  mrc_io_write_a3(io, mrc_a3_name(a3), a3);
}

static void
_mrc_a3_read(struct mrc_a3 *a3, struct mrc_io *io)
{
  a3->domain = (struct mrc_domain *)
    mrc_io_read_obj_ref(io, mrc_a3_name(a3), "domain", &mrc_class_mrc_domain);
  mrc_a3_setup(a3);
  mrc_io_read_a3(io, mrc_a3_name(a3), a3);
}

bool
mrc_a3_same_shape(struct mrc_a3 *a3_1, struct mrc_a3 *a3_2)
{
  if (a3_1->nr_comp != a3_2->nr_comp)
    return false;

  if (a3_1->sw != a3_2->sw)
    return false;

  int nr_patches_1, nr_patches_2;
  struct mrc_patch *patches_1 = mrc_domain_get_patches(a3_1->domain, &nr_patches_1);
  struct mrc_patch *patches_2 = mrc_domain_get_patches(a3_2->domain, &nr_patches_2);
  if (nr_patches_1 != nr_patches_2)
    return false;

  mrc_a3_foreach_patch(a3_1, p) {
    for (int d = 0; d < 3; d++) {
      if (patches_1[p].ldims[d] != patches_2[p].ldims[d])
	return false;
      if (patches_1[p].off[d] != patches_2[p].off[d])
	return false;
    }
  }
  return true;
}

// ----------------------------------------------------------------------
// mrc_class_mrc_a3

#define VAR(x) (void *)offsetof(struct mrc_a3, x)
static struct param mrc_a3_params_descr[] = {
  { "nr_comps"        , VAR(nr_comp)      , PARAM_INT(1)           },
  { "sw"              , VAR(sw)           , PARAM_INT(0)           },
  {},
};
#undef VAR

struct mrc_class_mrc_a3 mrc_class_mrc_a3 = {
  .name         = "mrc_a3",
  .size         = sizeof(struct mrc_a3),
  .param_descr  = mrc_a3_params_descr,
  .destroy      = _mrc_a3_destroy,
  .setup        = _mrc_a3_setup,
  .view         = _mrc_a3_view,
  .read         = _mrc_a3_read,
  .write        = _mrc_a3_write,
};

