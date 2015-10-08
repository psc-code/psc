#include <mrc_params.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrctest.h>

#include <assert.h>
#include <string.h>
#include <math.h>

const int sw = 2;

static struct mrc_fld *
setup_with_domain(struct mrc_domain *domain)
{
  struct mrc_fld *fld = mrc_domain_fld_create(domain, sw, "1,2,3,4,5,6");
  //struct mrc_fld *fld = mrc_domain_fld_create(domain, sw, NULL);
  //  mrc_fld_set_param_int_array(fld, "offs", 5, (int [5]) { 0, 2, 3, 0, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  //mrc_fld_view(fld);
  return(fld);
}

static struct mrc_fld *
setup_with_domain_nameless_comps(struct mrc_domain *domain)
{
  struct mrc_fld *fld = mrc_domain_fld_create(domain, sw, NULL);
  //  mrc_fld_set_param_int_array(fld, "offs", 5, (int [5]) { 0, 2, 3, 0, 0 });
  mrc_fld_set_from_options(fld);
  mrc_fld_set_param_int(fld, "nr_comps", 6);
  mrc_fld_setup(fld);
  //mrc_fld_view(fld);
  return(fld);
}

static void
init_standard(struct mrc_fld *fld)
{
  const char *fldtype = mrc_fld_type(fld);  
  int *offs = fld->_ghost_offs;
  int *dims = fld->_ghost_dims;
  for (int i4 = offs[4]; i4 < offs[4] + dims[4]; i4++) {
    for (int i3 = offs[3]; i3 < offs[3] + dims[3]; i3++) {
      for (int i2 = offs[2]; i2 < offs[2] + dims[2]; i2++) {
	for (int i1 = offs[1]; i1 < offs[1] + dims[1]; i1++) {
	  for (int i0 = offs[0]; i0 < offs[0] + dims[0]; i0++) {
	    if (strncmp(fldtype, "double", 6) == 0) {
	      MRC_D5(fld, i0,i1,i2,i3,i4) = i4 * 10000 + i3 * 1000 + i0 * 100 + i1 * 10 + i2;
	    }
	    else if (strncmp(fldtype, "float", 5) == 0) {
	      MRC_S5(fld, i0,i1,i2,i3,i4) = i4 * 10000 + i3 * 1000 + i0 * 100 + i1 * 10 + i2;
	    }
	    else if (strncmp(fldtype, "int", 3) == 0) {
	      MRC_FLD(fld, int, i0,i1,i2,i3,i4) = i4 * 10000 + i3 * 1000 + i0 * 100 + i1 * 10 + i2;
	    }
	  }
	}
      }
    }
  }
}


/* #define assert_equal(a, b, tol) { \ */
/*   double _norm = sqrt((a)*(a) - (b)*(b))/fabs(a); \ */
/*   if (_norm > (tol)) {				  \ */
/*     fprintf(stderr, "i0 %d i1 %d i2 %d i3 %d i4 %d f1 %g f2 %g: norm %g > %g\n", i0, i1, i2, i3, i4,(a),(b), _norm, tol); \ */
/*     assert(0); \ */
/*     } \ */
/*   } */

#define assert_equal(a, b, tol) assert((a) == (b))
static void
check_standard(struct mrc_fld *fld1, struct mrc_fld *fld2)
{
  const char *fldtype = mrc_fld_type(fld1); 
  int *offs = fld1->_offs.vals;
  int *dims = fld1->_dims.vals;
  for (int i4 = offs[4]; i4 < offs[4] + dims[4]; i4++) {
    for (int i3 = offs[3]; i3 < offs[3] + dims[3]; i3++) {
      for (int i2 = offs[2]; i2 < offs[2] + dims[2]; i2++) {
	for (int i1 = offs[1]; i1 < offs[1] + dims[1]; i1++) {
	  for (int i0 = offs[0]; i0 < offs[0] + dims[0]; i0++) {
	    if (strncmp(fldtype, "double",6) == 0) {
	      assert_equal(MRC_D5(fld1, i0,i1,i2,i3,i4), MRC_D5(fld2, i0,i1,i2,i3,i4), 1e-8);
	    }
	    else if (strncmp(fldtype, "float",5) == 0) {
	      assert_equal(MRC_S5(fld1, i0,i1,i2,i3,i4), MRC_S5(fld2, i0,i1,i2,i3,i4), 1e-5);
	    }
	    else if (strncmp(fldtype, "int",3) == 0) {
	      assert(MRC_FLD(fld1, int, i0,i1,i2,i3,i4) == MRC_FLD(fld2, int, i0,i1,i2,i3,i4));
	    }
	    else {
	      assert(0);
	    }
	  }
	}
      }
    }
  }

  assert(fld1->_nr_allocated_comp_name == fld2->_nr_allocated_comp_name);

  for (int m = 0; m < fld1->_nr_allocated_comp_name; m++) {
    const char *name1 = mrc_fld_comp_name(fld1, m);
    const char *name2 = mrc_fld_comp_name(fld2, m);
    // Check if both names are NULL (or the same memory location)
    if ( name1 != name2) {
      // Make sure neither one individually is NULL
      assert(name1);
      assert(name2);
      // Compare the actual text
      assert(strcmp(name1, name2) == 0);
    }
  }
}


int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "multi");
  mrc_domain_set_param_int3(domain, "m", (int[3]){12, 3, 4});
  mrc_domain_set_param_int3(domain, "np", (int[3]){1,1,1});
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  //  mrc_domain_view(domain);

  //struct mrc_fld *fld = setup_without_domain();
  bool nameless = false;
  mrc_params_get_option_bool("nameless", &nameless);

  struct mrc_fld *fld;
  if ( nameless ) {
    fld = setup_with_domain_nameless_comps(domain);
  } else {
    fld = setup_with_domain(domain);
  }

  init_standard(fld);


  struct mrc_io *io = mrc_io_create(mrc_fld_comm(fld));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  //mrc_io_view(io);


  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/fld", "test_fld", fld);
  mrc_io_close(io);

  mrc_io_destroy(io);

  io = mrc_io_create(mrc_fld_comm(fld));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  //mrc_io_view(io);

  mrc_io_open(io, "r", 0, 0.);
  struct mrc_fld *fld2 = mrc_io_read_path(io, "/fld", "test_fld", mrc_fld);
  mrc_io_close(io);
  mrc_io_destroy(io);

  // mrc_fld_view(fld2);

  assert(strcmp(mrc_fld_type(fld), mrc_fld_type(fld2)) == 0);
  const int *offs1 = mrc_fld_offs(fld);
  const int *offs2 = mrc_fld_offs(fld2);
  const int *dims1 = mrc_fld_dims(fld);
  const int *dims2 = mrc_fld_dims(fld2);
  const int *goffs1 = mrc_fld_ghost_offs(fld);
  const int *goffs2 = mrc_fld_ghost_offs(fld2);
  const int *gdims1 = mrc_fld_ghost_dims(fld);
  const int *gdims2 = mrc_fld_ghost_dims(fld2);
  for (int d = 0; d < 5; d++) {
    assert(offs1[d] == offs2[d]);
    assert(dims1[d] == dims2[d]);
    assert(goffs1[d] == goffs2[d]);
    assert(gdims1[d] == gdims2[d]);
  }
  
  check_standard(fld, fld2);

  mrc_fld_destroy(fld);
  mrc_fld_destroy(fld2);
  mrc_domain_destroy(domain);

  MPI_Finalize();
  return(0);
}
