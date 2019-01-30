
#include "psc_fields.h"
#include "psc_fields_as_single.h"
#include "fields.hxx"

using Fields = Fields3d<fields_t>;

#include "psc.h"

// ======================================================================
// helpers

// ----------------------------------------------------------------------
// make_mfields

static struct psc_mfields *
make_mfields(MPI_Comm comm)
{
  struct mrc_domain *domain = mrc_domain_create(comm);
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  
  struct psc_mfields *mflds = psc_mfields_create(comm);
  psc_mfields_set_type(mflds, FIELDS_TYPE);
  psc_mfields_set_param_obj(mflds, "domain", domain);
  psc_mfields_set_param_int(mflds, "nr_fields", 3);
  psc_mfields_set_param_int3(mflds, "ibn", (int [3]) { 1, 1, 1 });
  psc_mfields_set_from_options(mflds);
  psc_mfields_setup(mflds);

  // FIXME, I should be able to do this, but mflds still refers to it, and I'm
  // pretty sure there's no ref counting.
  //mrc_domain_destroy(domain);
  
  return mflds;
}

// ======================================================================
// the actual tests

// ----------------------------------------------------------------------
// test_create_destroy

static void
test_create_destroy(MPI_Comm comm)
{
  struct psc_mfields *mflds = psc_mfields_create(comm);
  psc_mfields_destroy(mflds);
}

// ----------------------------------------------------------------------
// test_setup

static void
test_setup(MPI_Comm comm)
{
  struct psc_mfields *mflds = make_mfields(comm);
  psc_mfields_view(mflds);
  psc_mfields_destroy(mflds);
}

// ----------------------------------------------------------------------
// test_print

static void
test_print(MPI_Comm comm)
{
  struct psc_mfields *mflds = make_mfields(comm);
  psc_mfields_view(mflds);

  mfields_t mf(mflds);
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = mf[p];
    Fields F(flds);

    for (int m = 0; m < flds.nr_comp; m++) {
      for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	  mprintf("mflds[%d, %d:%d,%d,%d, %d] =",
		  m, flds.ib[0], flds.ib[0] + flds.im[0], jy, jz, p);
	  for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	    printf(" %g", F(m, jx,jy,jz));
	  }
	  printf("\n");
	}
	printf("\n");
      }
      printf("\n");
    }
  }
  
  psc_mfields_destroy(mflds);
}

// ----------------------------------------------------------------------
// dummy_psc_main
//
// FIXME, just a hack

void
dummy_psc_main()
{
  // FIXME: If we don't reference psc_main, we get link errors
  // from missing global variables (for timing etc),
  psc_main(NULL, NULL, NULL);
}

// ----------------------------------------------------------------------
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  mprintf("-- test_create_destroy\n");
  test_create_destroy(comm);

  mprintf("-- test_setup\n");
  test_setup(comm);
  
  mprintf("-- test_print\n");
  test_print(comm);
  
  MPI_Finalize();
}
