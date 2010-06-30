#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psc.h"
#include "psc_cbe.h"
#include "psc_alf.h"
#include "util/profile.h"

// The following lines are not my own. They are a standard in the IBM examples, and who
// am I to think I can do better? --steve

/// Max length of path to PPU image
///
/// Not needed for Cell only systems, but necessary for hybrid.
#define PATH_BUF_SIZE 1024 


/// Name of spu shared object library
char library_name[PATH_BUF_SIZE];


/// Complete path to SPU image
char spu_image_path[PATH_BUF_SIZE];   

/// Name of binary spu image embeded in so library
const char *spu_image_name = "alf_cbe_push_part_yz_kernel"; 

/// Name of computational kernel function
///
/// Specified by ALF_EXPORT_API in computation kernel file
const char *kernel_name = "cbe_push_part_yz_kernel";

/// Name of accelerator side DTL prep function.
///
/// Specified by ALF_EXPORT_API in computation kernel file
const char *prep_in_name = "push_yz_prep_in_dtl";

/// Name of accelerator side DTL out prep function.
///
/// Specified by ALF_EXPORT_API in computation kernel file
const char *prep_out_name = "push_yz_prep_out_dtl";

// If we want to use this on hybrid system (ie Roadrunner) we would need a variable for the 
// PPE image path.

/// Data common to all task instances
///
/// Contains global simulation parameters and coeff
push_task_context_t tsk_ctx __attribute__((aligned(128))); 


////////////////////////////
/// Cell optimized implementation of the YZ particle pusher.
///
/// This function serves as the entry point to the particle pusher, and runs on the PPE.
/// Currently it initializes the Accelerated Library Framework (ALF) environment and 
/// initializes the task context. We really should transition to setting up a shared 
/// ALF environment at the start of a run, and only initializing the task context once 
/// in the course of the run, as this process is fairly expensive. It might be possible
/// to speed it up by directly editing the ALF variables instead of using IBM's function
/// calls, but I really think it's better if we only have to do this once per run.
///
/// After the task is initialized, the PPE does a domain decomposition based on the grid,
/// currently only one grid cell at a time. It then passes this information to the SPUs, which
/// each perform their own domain decomposition on the particles to account for the limitations of
/// the SPUs. ALF handles all the data transfer and work block scheduling. When the SPUs are 
/// finished running, control is returned to the PPE and the pusher exits. We could probably find
/// something productive for the PPE to be doing during this time...
static void 
do_cbe_push_part_yz(particles_cbe_t *pp, fields_cbe_t *pf){
  
  int ierr;
  alf_task_desc_handle_t push_desc_handle;
  alf_task_handle_t push_handle;
  alf_wb_handle_t wb_push_handle;
  push_wb_context_t wb_ctx;

  // Max number of particles on the spu at one time.
  unsigned int maxpart = 100000/(int)sizeof(particle_cbe_t); 

  sprintf(library_name, "cbe_psc_spu.so");
  
  ///< \FIXME Max part per WB should probably be calculated at runtime. There is a serious problem at the moment with some sort of memory fault on the spus that limits the max number of particles to 200. I think I have a handle on what's going on, and will hopefully fix it soon. 

  // I think the 0 here means spe...
  ierr = alf_task_desc_create(alf_env_shared, 0, &push_desc_handle); ACE;
  
  // Don't need to set type, as default is WB

  // Tell ALF we're partitioning on accelerator
  ierr = alf_task_desc_set_int32(push_desc_handle, ALF_TASK_DESC_PARTITION_ON_ACCEL,
				 1); ACE;

  //Set size of task context buffer
  ierr = alf_task_desc_set_int32(push_desc_handle, ALF_TASK_DESC_TSK_CTX_SIZE,
				 sizeof(push_task_context_t)); ACE;
  
  //Set size of param buffer for each WB
  ierr = alf_task_desc_set_int32(push_desc_handle, ALF_TASK_DESC_WB_PARM_CTX_BUF_SIZE, 
				 sizeof(push_wb_context_t) +
				 (2 * psc.ibn[1] + 1) * (2*psc.ibn[2] + 1)
				 * 4 * sizeof(fields_cbe_real_t) 
				 + (16-sizeof(push_wb_context_t)%16));ACE;// Number of cells + ghosts); ACE;

  //Set size of input buffer (will hold E/B fields)
  ierr= alf_task_desc_set_int32(push_desc_handle, ALF_TASK_DESC_WB_IN_BUF_SIZE,
  				(2 * psc.ibn[1] + 1) * (2*psc.ibn[2] + 1) // Number of cells + ghosts
				* 6 * sizeof(fields_cbe_real_t)); ACE;
			       
  //Set size of output buffer (will hold current densities)
  ierr = alf_task_desc_set_int32(push_desc_handle, ALF_TASK_DESC_WB_OUT_BUF_SIZE,
  				(2 * psc.ibn[1] + 1) * (2*psc.ibn[2] + 1) // Number of cells + ghosts
				* 4 * sizeof(fields_cbe_real_t)); ACE;

  //Set size of in/out buffer (will hold particles, should be set to max allowable particles/WB)
  ierr = alf_task_desc_set_int32(push_desc_handle, ALF_TASK_DESC_WB_INOUT_BUF_SIZE, 
				 (maxpart+1) * sizeof(particle_cbe_t)); ACE;

  ierr = alf_task_desc_set_int32(push_desc_handle, ALF_TASK_DESC_NUM_DTL, 
				 20); ACE;

  ierr = alf_task_desc_set_int32(push_desc_handle, ALF_TASK_DESC_NUM_DTL_ENTRIES, 
				 200); ACE;

  //There are a few more things that could be set here pertaining to DTLs and accelerator stacks.
  // At the moment, I don't think it's neccessary to change them from the default, but if it 
  // does become needed, this is where those options belong.


  // Set the name of the library containing the image of the computational kernel
  ierr = alf_task_desc_set_int64(push_desc_handle, ALF_TASK_DESC_ACCEL_LIBRARY_REF_L,
				 (unsigned long long) library_name); ACE;

  // Set the name of the accelerator image
  ierr = alf_task_desc_set_int64(push_desc_handle, ALF_TASK_DESC_ACCEL_IMAGE_REF_L, 
				 (unsigned long long) spu_image_name); ACE;

  // Set the name of the computational kernel function
  ierr = alf_task_desc_set_int64(push_desc_handle, ALF_TASK_DESC_ACCEL_KERNEL_REF_L, 
				 (unsigned long long) kernel_name); ACE;

  // Set the name of the accelerator dtl-in prepare function
  ierr = alf_task_desc_set_int64(push_desc_handle, ALF_TASK_DESC_ACCEL_INPUT_DTL_REF_L,
				 (unsigned long long) prep_in_name); ACE;

  // Set the name of the accelerator dtl-out prepare function
  ierr = alf_task_desc_set_int64(push_desc_handle, ALF_TASK_DESC_ACCEL_OUTPUT_DTL_REF_L,
				 (unsigned long long) prep_out_name); ACE;

  // Init params and register in handle
  psc_alf_init_push_task_context(&tsk_ctx, &push_desc_handle, pf);

  // Done with description, time to create the task.
  ierr = alf_task_create(push_desc_handle, &tsk_ctx, 
			 0, // Use all available spes
			 0, // Let alf handle the scheduling
			 0, // No bundling of work blocks at this moment
			 &push_handle); ACE;


  // Creating a dataset will neither hurt nor help on a straight Cell system, but
  // on a hybrid system it will improve scheduling of data transfer between the host
  // and the accelerator. Might as well do it now.
#if 0
  alf_dataset_handle_t dataset_h;
  
  ierr = alf_dataset_create(push_handle, &dataset_h); ACE;
  ierr = alf_dataset_buffer_add(dataset_h, pp->particles, 
				psc.pp.n_part * sizeof(particle_cbe_t),
				ALF_DATASET_READ_WRITE);
  ierr = alf_task_dataset_associate(push_handle, dataset_h); ACE;
#endif

  int n = 0; // Important! This assumes particles are sorted!!! Otherwise, life will be bad!
             // Okay, not really bad. You'll just end up doing more workblocks than needed, 
             // and the code will run slow. 

  int num_wb = 0; // count the number of workblocks for debugging purposes. 
  // This is ugly, and there must be a better way. 
  // Splitting into wb bigger than one cell is a good idea. 
  int avgncell = psc.pp.n_part / ((psc.ihi[2]-psc.ilo[2]) * (psc.ihi[1] - psc.ilo[1]));
  int nend;
  
  // Need a pointer that moves through the padding array. For a variety of
  // reasons, very bad things happen if you try to use a single memory location
  // to pad out the particles. I'll write more about it later.
  particle_cbe_t * padding = pp->null_particles;

  spu_curr_cache_t * current_cache;
  current_cache = calloc(psc.fld_size, sizeof(spu_curr_cache_t));
  
  while( n < psc.pp.n_part){

    nend = n + avgncell - 1;
    
    cbe_real cni_ref = pp->particles[n].cni;
    
    if(nend > (psc.pp.n_part -1)) nend = psc.pp.n_part - 1;
    
    cbe_real cni_curr = pp->particles[nend].cni;
    
    if(fabs(cni_ref - cni_curr) < 0.5) {
      while(fabs(cni_ref - cni_curr) < 0.5) {
	nend++;
	if(nend >= psc.pp.n_part) break;
	cni_curr = pp->particles[nend].cni;
      }
      nend--;
    }
    else {
      while(fabs(cni_ref - cni_curr) >= 0.5){
	nend--;
	cni_curr = pp->particles[nend].cni;
      }
    }
    // FIXME we need to add some empty particles on to the end of the in/out buffer to pad 
    // out the vector ops. The 10 above is necessary for a poorly understood SPU problem.

    // Need to do fields here

    wb_ctx.p_start = (unsigned long long) &(pp->particles[n]);    
    wb_ctx.p_pad = (unsigned long long) padding;
    // Assign values in WB for padding
    cbe_real yref = floor(pp->particles[n].yi), zref = floor(pp->particles[n].zi);
    for(int i = 0; i < VEC_SIZE; i++){
      padding->xi = pp->particles[n].xi;
      padding->yi = pp->particles[n].yi;
      padding->zi = pp->particles[n].zi;
      padding->pxi = 0.;
      padding->pyi = 0.;
      padding->pzi = 0.;
      padding->qni = 1.;
      padding->mni = 1.;
      padding->wni = 0.;
      padding->cni = 0.;
      padding++;
    }
    wb_ctx.lo_npart = nend - n + 1;
    wb_ctx.maxpart = maxpart;
    // NB these aren't really lg for x, they're lo instead. 
    // I don't really care about the ghost points in x.
    wb_ctx.wb_lg[0] = psc.ilo[0];
    wb_ctx.wb_hg[0] = psc.ihi[0]-1;
    wb_ctx.wb_lg[1] = (int)(yref/psc.dx[1]) - psc.ibn[1];
    wb_ctx.wb_hg[1] = (int)(yref/psc.dx[1]) + psc.ibn[1]; 
    wb_ctx.wb_lg[2] = (int)(zref/psc.dx[2]) - psc.ibn[2];
    wb_ctx.wb_hg[2] = (int)(zref/psc.dx[2]) + psc.ibn[2];


    int iterations;
    if( (wb_ctx.lo_npart % maxpart) == 0) {
      iterations = wb_ctx.lo_npart / maxpart;
    } else {
      iterations = wb_ctx.lo_npart / maxpart + 1;
    }

    wb_current_cache_init(&current_cache[num_wb], &wb_ctx);

    ierr = alf_wb_create(push_handle, ALF_WB_MULTI, iterations, &wb_push_handle); ACE;
    
    alf_wb_parm_add(wb_push_handle, (void *)(&wb_ctx), sizeof(push_wb_context_t), ALF_DATA_BYTE, 0);
    
    ierr = alf_wb_enqueue(wb_push_handle); ACE;
    
    n = nend + 1;
    num_wb ++;
  }

  //  printf("Number of WB: %d\n", num_wb);
  // Starts automatically once finalized
  ierr = alf_task_finalize(push_handle); ACE;
  
  // For now, the PPE just waits
  ierr = alf_task_wait(push_handle, -1); ACE;
  
  ierr = alf_task_destroy(push_handle); ACE;
  for(int m = 0; m < num_wb; m++){ 
    wb_current_cache_store(pf, &current_cache[m]);
    fflush(stdout);
    free(current_cache[m].flds);
    current_cache[m].flds = NULL;
  }
  free(current_cache);
  current_cache = NULL;
  
}

void cbe_push_part_yz()
{
  particles_cbe_t pp;
  fields_cbe_t pf;
  particles_cbe_get(&pp);
  fields_cbe_get(&pf, EX, EX + 6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("cbe_part_yz", 1., 0, psc.pp.n_part * 9 * sizeof(cbe_real));
  }
  prof_start(pr);
  do_cbe_push_part_yz(&pp, &pf);
  prof_stop(pr);

  particles_cbe_put(&pp);
  fields_cbe_put(&pf, JXI, JXI + 3);
}
