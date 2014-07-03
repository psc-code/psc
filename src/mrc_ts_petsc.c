
#include <mrc_ts_private.h>
#include <math.h>
#include <petscts.h>

#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_ddc.h>

#include <assert.h>

#define CE assert(ierr == 0)

#define pfb PetscFunctionBegin
#define pfr PetscFunctionReturn(0)

#define mrc_ts_petsc(x) mrc_to_subobj(x, struct mrc_ts_petsc)


// Some typedefs to make the method fetching a little easier
typedef Vec (*fgp_t)(struct mrc_fld *);
typedef void (*fpp_t)(struct mrc_fld *, Vec *);
typedef void (*fsp_t)(struct mrc_fld *, Vec);

static int 
wrap_mrc_monitors(TS petsc_ts, int nstep, PetscReal t, Vec xg, void *_ctx)
{
  pfb;
  // I didn't really think about this all that well with the get/set petsc_vec
  // bits. In essence, this doesn't do anything except change the gotten_petsc_vec
  // tag. I suppose the time stepper is sort of a special case...
  struct mrc_ts *ts = (struct mrc_ts *) _ctx;
  struct mrc_ts_petsc *sub = mrc_ts_petsc(ts);
  fpp_t fld_put_petsc = (fpp_t) mrc_fld_get_method(sub->xg, "put_petsc_vec");
  fld_put_petsc(sub->xg, &xg);
  struct mrc_ddc *ddc = mrc_domain_get_ddc(sub->xg->_domain);
  mrc_ddc_global_to_local_fld(ddc, sub->xg, (struct mrc_fld *) ts->x);
  // This needs to get updated here, rather than in the rhs wrapper
  ts->n = nstep;
  
  mrc_ts_monitors_run(ts);
  // Let's set a variable that is about to go out of scope. Yah for uselessness!
  fgp_t fld_get_petsc = (fgp_t) mrc_fld_get_method(sub->xg, "get_petsc_vec");
  xg = fld_get_petsc(sub->xg);
  pfr;
}

// this is almost a duplicate except it changes the number of ghosts.
// FIXME: Adapt this to work without domains
static inline void
set_global_from_fld(struct mrc_fld *gfld, struct mrc_fld *from)
{
  mrc_fld_set_type(gfld, mrc_fld_type(from));
  mrc_fld_set_param_obj(gfld, "domain", from->_domain);
  mrc_fld_set_param_int_array(gfld, "offs",from->_offs.nr_vals, from->_offs.vals);
  mrc_fld_set_param_int(gfld, "nr_comps", mrc_fld_nr_comps(from));
  mrc_fld_set_param_int(gfld, "nr_spatial_dims", from->_nr_spatial_dims);
  mrc_fld_set_param_int(gfld, "nr_ghosts", 0);
  mrc_fld_set_param_int_array(gfld, "sw", from->_sw.nr_vals, NULL);
}

static int 
wrap_mrc_rhs_function(TS petsc_ts, PetscReal t, Vec xg, Vec F, void *_ctx)
{
  pfb;
  struct mrc_ts *ts = (struct mrc_ts *) _ctx;

  struct mrc_ts_petsc *sub = mrc_ts_petsc(ts);

  // I don't really like it, but we don't know if the vectors coming in are
  // ts->x and sub->F or intermediate values, which is really a pain
  // So we need to create mrc_fld wrappers every time, which means
  // we are aliasing the petsc vectors. Blah.
  struct mrc_fld *F_fld = mrc_fld_create(mrc_fld_comm(sub->F));
  set_global_from_fld(F_fld, sub->F);
  fsp_t fld_set_petsc = (fsp_t) mrc_fld_get_method(sub->F, "set_petsc_vec");
  fld_set_petsc(F_fld, F);  
  mrc_fld_setup(F_fld);

  struct mrc_fld *xg_fld = mrc_fld_create(mrc_fld_comm(sub->xg));
  set_global_from_fld(xg_fld, sub->xg);
  fld_set_petsc(xg_fld, xg);  
  mrc_fld_setup(xg_fld);
  
  ts->time = t;

  // we need to get the info from the global state vector xg (which petsc wants) back into the local
  // state vector x (which libmrc wants). Fun times.
  
  struct mrc_ddc *ddc = mrc_domain_get_ddc(sub->xg->_domain);
  mrc_ddc_global_to_local_fld(ddc, xg_fld, sub->xlocal);
  
  ts->rhsf(ts->rhsf_ctx, 
	   mrc_fld_to_mrc_obj(F_fld), 
	   (float) t, 
	   mrc_fld_to_mrc_obj(sub->xlocal));
  
  mrc_fld_destroy(xg_fld);
  mrc_fld_destroy(F_fld);
  pfr;
}

// wrap_mrc_jacobian_function
//
// Wrap Jacobian calc. Add in some extra bits to only recompute it when
// we really need to (to improve efficiency)
// Right now I can't think of anyway to write a jacobian function
// without some matrix format, so I'll just have to lock calc_jac 
// into petsc_enabled stuff for now.

// Set both the calc_jac function AND the get_jac, which will only be called at setup time.
void
mrc_ts_petsc_set_jacobian_function(struct mrc_ts *ts,
			     int (*calc_jac)(Mat J, Vec x, float t, void *ctx),
			     int (*get_jac_mat)(void *ctx, Mat *M),
			     int (*calc_pre_mat)(Mat Pre, Vec x, float t, void *ctx),
			     int (*get_pre_mat)(void *ctx, Mat *M),
			     void *ctx)
{
  struct mrc_ts_petsc *sub = mrc_ts_petsc(ts);
  sub->calc_jac = calc_jac;
  sub->get_jac_mat = get_jac_mat;
  sub->calc_pre_mat = calc_pre_mat;
  sub->get_pre_mat = get_pre_mat;
  sub->jacf_ctx = ctx;
}


static int
wrap_mrc_jacobian_function(TS ts, PetscReal t, Vec xg, Mat *J, Mat *B,
			   MatStructure *flag, void *ctx)
{
  struct mrc_ts *tl = ctx;
  struct mrc_ts_petsc *sub = mrc_ts_petsc(tl);
  int ierr, step, it, freq, freq_step;
  SNES snes;
  int bsin, bsctx, sizein, sizectx;
  ierr = VecGetBlockSize(xg, &bsin); CE;
  ierr = VecGetLocalSize(xg, &sizein); CE;
  ierr = VecGetBlockSize(sub->xg_vec, &bsctx); CE;
  ierr = VecGetLocalSize(sub->xg_vec, &sizectx); CE;
  if ( sizein == sizectx && bsin != bsctx) {
    VecSetBlockSize(xg, bsctx);
  }

  PetscFunctionBegin;
  ierr = TSGetSNES(ts, &snes); CE;
  ierr = SNESGetIterationNumber(snes, &it); CE;
  ierr = TSGetTimeStepNumber(ts, &step); CE;
  // !!!
  freq = 7;
  // FIXME at the very least, make that from the last recomputation...
  freq_step = 20000;

  /* static Mat Bsave; */
  /* if (!Bsave) { */
  /*   ierr = MatDuplicate(*J, MAT_DO_NOT_COPY_VALUES, &Bsave); CE; */
  /* } */

  if ((step % freq_step) == 0 && it == 0) {
    goto compute;
  }
  // if the snes takes more than 7 iterations to converge, recompute
  // the jacobian.
  if (((it % freq) == 0 && it != 0)) {
    goto compute;
  }


  // Modifications Kai made to this timestepper based on changes in Petsc.
  // The petsc timestepper might be in need of some cleaning sometime soon...
  goto compute;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Skipping Jacobian...\n"); CE;
  /* assert(Bsave); */
  /* ierr = MatCopy(Bsave, *B, SAME_NONZERO_PATTERN); CE; */
#warning should not be necessary to copy jacobian back???
  *flag = SAME_PRECONDITIONER;
  goto out;

 compute:
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Computing Jacobian... "); CE;
    //    ierr = LogEventBegin(KG_calc_jac); CE;
    ierr = sub->calc_jac(*J, xg, (float) t, sub->jacf_ctx); CE;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CE;
    //    ierr = LogEventEnd(KG_calc_jac); CE;
    *flag = SAME_NONZERO_PATTERN;
    if (sub->sep_pre) {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Computing Preconditioner Matrix... "); CE;
      ierr = sub->calc_pre_mat(*B, xg, (float) t, sub->jacf_ctx); CE;
      *flag = *((MatStructure *) sub->pre_mat_structure);
      ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CE;
    }
    //    ierr = MatView(*J, PETSC_VIEWER_STDOUT_WORLD); CE;
    //    ierr = print_jac(*J); CE;
    //ierr = MatCopy(*B, Bsave, SAME_NONZERO_PATTERN); CE;
#if 0
    {
      int i,j, jx,jy,m, ix = 0, iy = 4, im = B0;
      double val, val2;

      i = (iy*128+ix)*8 + im;
      for (jx = ix-1 >= 0?:0; jx <= ix+1; jx++) {
	for (jy = iy-1; jy <= iy+1; jy++) {
	  for (m = 0; m < 8; m++) {
	    j = (jy*128+jx)*8 + m;
	    ierr = test_jac(tl, *J, xg, i, j, &val, &val2);
	    if (val || val2) {
	      printf("(%d,%d,%d) (%d,%d,%d) %20g %20g %s\n", ix,iy,im, jx,jy,m,
		     val, val2, fabs(val-val2)/MAX(fabs(val),fabs(val2)) > 1e-7 ? "!!!" : "");
	    }
	  }
	}
	for (jy = iy-1+64; jy <= iy+1+64; jy++) {
	  for (m = 0; m < 8; m++) {
	    j = (jy*128+jx)*8 + m;
	    ierr = test_jac(tl, *J, xg, i, j, &val, &val2);
	    if (val || val2) {
	      printf("(%d,%d,%d) (%d,%d,%d) %20g %20g %s\n", ix,iy,im, jx,jy,m,
		     val, val2, fabs(val-val2) > 1e-10 ? "!!!" : "");
	    }
	  }
	}
      }
    }
#endif

 out:
  PetscFunctionReturn(0);
}


static void 
_mrc_ts_petsc_create(struct mrc_ts *ts)
{
  struct mrc_ts_petsc *sub = mrc_ts_petsc(ts);
  int ierr = TSCreate(PETSC_COMM_WORLD, &sub->petsc_ts); CE;
}

static void 
_mrc_ts_petsc_setup(struct mrc_ts *ts)
{
  struct mrc_ts_petsc *sub = mrc_ts_petsc(ts);
  assert(ts->x);
  struct mrc_fld *x_fld = (struct mrc_fld *) ts->x;
  int ierr;
  
  ierr = TSSetProblemType(sub->petsc_ts, TS_NONLINEAR); CE;

  // Create a global version of the state vector
  sub->xg = mrc_fld_create(mrc_fld_comm(x_fld));
  set_global_from_fld(sub->xg, x_fld);
  mrc_fld_setup(sub->xg);
  
  // Create the F vector. Should be almost identical to x, except
  // it will be global (no matter what), so we can't just use duplicate
  sub->F = mrc_fld_create(mrc_fld_comm(x_fld));
  set_global_from_fld(sub->F, x_fld);
  mrc_fld_setup(sub->F);

  fgp_t fld_get_petsc = (fgp_t) mrc_fld_get_method(sub->F, "get_petsc_vec");
  sub->F_vec = fld_get_petsc(sub->F);

  sub->xlocal = (struct mrc_fld *) ts->vec_duplicate(ts->x);

  sub->xg_vec = fld_get_petsc(sub->xg);
  // Set the global state vector
  ierr = TSSetSolution(sub->petsc_ts, sub->xg_vec); CE;

  // Register the RHS functions, and the vector to hold it.
  ierr = TSSetRHSFunction(sub->petsc_ts, sub->F_vec, wrap_mrc_rhs_function, (void *) ts); CE;

  // Register the monitor functions.
  ierr = TSMonitorSet(sub->petsc_ts, wrap_mrc_monitors, (void *) ts, PETSC_NULL); CE;



  // if we have a Jacobian, register that too.
  if (sub->calc_jac){
    assert(sub->get_jac_mat);
    sub->get_jac_mat(sub->jacf_ctx, &sub->J);
    if (sub->sep_pre) {
      assert(sub->get_pre_mat);
      assert(sub->calc_pre_mat);
      if (!sub->pre_mat_structure) {
	sub->pre_mat_structure = calloc(1, sizeof(MatStructure));
      }
      sub->get_pre_mat(sub->jacf_ctx, &sub->Pre);
      ierr = TSSetRHSJacobian(sub->petsc_ts, sub->J, sub->Pre,
			      wrap_mrc_jacobian_function, ts); CE;
    }
    else {
      ierr = TSSetRHSJacobian(sub->petsc_ts, sub->J, sub->J,
			      wrap_mrc_jacobian_function, ts); CE;
    }
  }


  ierr = TSSetTime(sub->petsc_ts, ts->time); CE;

  // See if there are any petsc specific options that the timestepper
  // needs to know about.
  ierr = TSSetFromOptions(sub->petsc_ts); CE;

  
}

static void 
_mrc_ts_petsc_destroy(struct mrc_ts *ts)
{
  struct mrc_ts_petsc *sub = mrc_ts_petsc(ts);
  // clean up all my helper fields
  mrc_fld_destroy(sub->F);
  mrc_fld_destroy(sub->xg);
  mrc_fld_destroy(sub->xlocal);
  int ierr;
  if (sub->J) {
    ierr = MatDestroy(&(sub->J));
  }
  ierr = TSDestroy(&(sub->petsc_ts)); CE;  
  
  if (sub->pre_mat_structure) free(sub->pre_mat_structure);

}

static void 
_mrc_ts_petsc_solve(struct mrc_ts *ts)
{
  struct mrc_ts_petsc *sub = mrc_ts_petsc(ts);

  struct mrc_fld *x = (struct mrc_fld *) ts->x;

  // copy the local vector x into the global vector xg before starting the petsc solver.
  // FIXME: HACKY and only works for 5D flds!
  const int *dims = mrc_fld_dims(x);
  for (int i4 = 0; i4 < dims[4]; i4++) {
    for (int i3 = 0; i3 < dims[3]; i3++) {
      for (int i2 = 0; i2 < dims[2]; i2++) {
	for (int i1 = 0; i1 < dims[1]; i1++) {
	  for (int i0 = 0; i0 < dims[0]; i0++) {
	    if (x->_data_type == MRC_NT_DOUBLE) {
	      MRC_D5(sub->xg, i0, i1, i2, i3, i4) = MRC_D5(x, i0, i1, i2, i3, i4);
	    } else if (x->_data_type == MRC_NT_FLOAT) {
	      MRC_S5(sub->xg, i0, i1, i2, i3, i4) = MRC_S5(x, i0, i1, i2, i3, i4);
	    } else {
	      assert(0);
	    }
	  }
	}
      }
    }
  }

  // Run solver here
  int ierr = TSSolve(sub->petsc_ts, NULL); CE;

  fpp_t fld_put_petsc = (fpp_t) mrc_fld_get_method(sub->xg, "put_petsc_vec");
  fld_put_petsc(sub->xg, &sub->xg_vec);
  fld_put_petsc(sub->F, &sub->F_vec);

  // now that solver is done, do a final copy of global into local
  struct mrc_ddc *ddc = mrc_domain_get_ddc(sub->xg->_domain);
  mrc_ddc_global_to_local_fld(ddc, sub->xg, (struct mrc_fld *) (ts->x));  

}


static struct mrc_obj_method mrc_ts_petsc_methods[] = {
  MRC_OBJ_METHOD("set_jac", mrc_ts_petsc_set_jacobian_function),
  {}
};

#define VAR(x) (void *)offsetof(struct mrc_ts_petsc, x)
static struct param mrc_ts_petsc_param_descr[] = {
  { "sep_pre"           , VAR(sep_pre)           , PARAM_BOOL(false),
  .help = "the preconditioner will be calculated from a matrix other than the jacobian" },
  { "pre_mat_structure" , VAR(pre_mat_structure) , PARAM_PTR(NULL),
    .help = "A pointer to the petsc structure flag describing the preconditioner matrix. NOT FOR COMMAND LINE USE!" },
  {},
};
#undef VAR


struct mrc_ts_ops mrc_ts_petsc_ops = {
  .name             = "petsc",
  .size             = sizeof(struct mrc_ts_petsc),
  .param_descr      = mrc_ts_petsc_param_descr,
  .create           = _mrc_ts_petsc_create,
  .setup            = _mrc_ts_petsc_setup,
  .destroy          = _mrc_ts_petsc_destroy,
  .solve            = _mrc_ts_petsc_solve,
  .methods          = mrc_ts_petsc_methods,
};
