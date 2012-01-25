/*
 *  psc_output_particles_c.c
 *  psc
 *
 *  Created by Nils Moschuering on 26.09.11.
 *  Copyright 2011 LMU. All rights reserved.
 *
 *  Hacked by Simon Jagoda on 11.11.11.
 */



#include "psc_output_particles_private.h"
#include "psc_output_particles_custom_binary.h"

#include "mrc_io.h"
#include "mpi.h"
#include "psc.h"
#include "mrc_domain.h"
#include "mrc_fld.h"
#include <mrc_params.h>
#include <assert.h>

#include "hdf5.h"

#include <mrc_profile.h>
#include <string.h>

#define to_psc_output_particles_custom_binary(out) ((struct psc_output_particles_custom_binary *)((out)->obj.subctx))

/**
	
 * Ansi C "itoa" based on Kernighan & Ritchie's "Ansi C":
	
 */
	
void strreverse(char* begin, char* end) {
	char aux;
	while(end>begin)
		aux=*end, *end--=*begin, *begin++=aux;}
	
void itoa(int value, char* str, int base) {
	static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
	char* wstr=str;
	int sign;
	if (base<2 || base>35){ *wstr='\0'; return; }
	if ((sign=value) < 0) value = -value;
	do *wstr++ = num[value%base]; while(value/=base);
	if(sign<0) *wstr++='-';
	*wstr='\0';
	strreverse(str,wstr-1);}


// ----------------------------------------------------------------------
// psc_output_particles_c_run

static void
psc_output_particles_custom_binary_run(struct psc_output_particles *out,
                              mparticles_base_t *particles)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_out_part", 1., 0, 0);
  }
  prof_start(pr);
  
  struct psc *psc = ppsc;
  struct psc_output_particles_custom_binary *out_c = to_psc_output_particles_custom_binary(out);
  if (psc->timestep == 0) {create_output_file(out_c);}
  if (psc->timestep >= out_c->next) {
    out_c->next += out_c->step;
    


    write_particles_to_file(out_c,particles);	
    
 
  }
  
  prof_stop(pr);
}


static void create_output_file(struct psc_output_particles_custom_binary *out_c) {
  struct psc *psc = ppsc;
  int node;
  MPI_Comm_rank(MPI_COMM_WORLD, &node);
  char filename[80];
  strcpy(filename, out_c->data_dir);//directory...
  strcat(filename,"/part"); 
  char buffer[5];
  itoa(node,buffer,10);
  strcat(filename,buffer);	//node...
  FILE *file = fopen(filename,"w");	//create output file
  
  int nodes=psc->domain.np[0]+psc->domain.np[1]*psc->domain.np[2];	//file header
  fwrite(&nodes,sizeof(int),1,file);
  fwrite(&psc->prm.nmax,sizeof(int),1,file);
  fwrite(&out_c->write_x,sizeof(bool),1,file);
  fwrite(&out_c->write_y,sizeof(bool),1,file);
  fwrite(&out_c->write_z,sizeof(bool),1,file);
  fwrite(&out_c->write_px,sizeof(bool),1,file);
  fwrite(&out_c->write_py,sizeof(bool),1,file);
  fwrite(&out_c->write_pz,sizeof(bool),1,file);
  fwrite(&out_c->write_q,sizeof(bool),1,file);
  fwrite(&out_c->write_m,sizeof(bool),1,file);
  
  fclose(file);


}


static void write_particles_to_file(struct psc_output_particles_custom_binary *out_c, mparticles_base_t *particles_in)	//actually particles to filesystem now
{  
  mparticles_t *particles = psc_mparticles_get_cf(particles_in, 0);
  struct psc *psc = ppsc;
  
  //generate filename
  int node;
  MPI_Comm_rank(MPI_COMM_WORLD, &node);
  char filename[80];
  strcpy(filename, out_c->data_dir);//directory...
  strcat(filename,"/part"); 
  char buffer[5];
  itoa(node,buffer,10);
  strcat(filename,buffer);	//node...
  FILE *file = fopen(filename,"a");	//finally open the damn thing
  
 
  
  fwrite(&psc->timestep,sizeof(int),1,file);	//current timestep
  
  for (int p = 0; p < particles->nr_patches; p++) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    
    //number of particles on this patch
    if(out_c->filter_func) {		//if filtering is active
    	int npart=0;
    	for(int i=0; i<pp->n_part; i++) {
    		if(out_c->filter_func(particles_get_one(pp, i))==true) {npart++;}
    	}
    fwrite(&npart,sizeof(int),1,file);
    }
    
    else{fwrite(&pp->n_part,sizeof(int),1,file);}
    
    //------------------------------------------------------------------
    //testing start



    //testing end
    //------------------------------------------------------------------
    
    for (int n = 0; n < pp->n_part; n++) {
      particle_t *part = particles_get_one(pp, n);
      //write particle data
      if(out_c->filter_func && !(out_c->filter_func(part))) continue;
      if(out_c->write_x==true){fwrite(&part->xi,sizeof(double),1,file);}
      if(out_c->write_y==true){fwrite(&part->yi,sizeof(double),1,file);}
      if(out_c->write_z==true){fwrite(&part->zi,sizeof(double),1,file);}
      if(out_c->write_px==true){fwrite(&part->pxi,sizeof(double),1,file);}
      if(out_c->write_py==true){fwrite(&part->pyi,sizeof(double),1,file);}
      if(out_c->write_pz==true){fwrite(&part->pzi,sizeof(double),1,file);}
      if(out_c->write_q==true){fwrite(&part->qni,sizeof(double),1,file);}
      if(out_c->write_m==true){fwrite(&part->mni,sizeof(double),1,file);} 
     
      
      
 //if(out_c->filter_func==true) printf("true \n");
  //if(out_c->filter_func==false) printf("false \n");
    }
  }
  fclose(file);
}

void
psc_output_particles_custom_binary_setfilter(struct psc_output_particles *out, bool (*new_filter_func)(particle_t *part))
{
  struct psc_output_particles_custom_binary *out_c = to_psc_output_particles_custom_binary(out);
  out_c->filter_func=new_filter_func;
}

static void
psc_output_particles_custom_binary_create(struct psc_output_particles *out)
{
  struct psc_output_particles_custom_binary *out_c = to_psc_output_particles_custom_binary(out);
  out_c->io = mrc_io_create(psc_output_particles_comm(out));
  //I'd like to have the outdir from mrc_io as the default value of data_dir
  mrc_io_set_from_options(out_c->io);
  const char *outdir=NULL;
  mrc_obj_get_param_string((struct mrc_obj*)out_c->io, "outdir" ,&outdir);
  psc_output_particles_set_param_string(out, "data_dir" ,outdir);
}

static void
psc_output_particles_custom_binary_setup(struct psc_output_particles *out)
{
  struct psc_output_particles_custom_binary *out_c = to_psc_output_particles_custom_binary(out);

  const char *outdir=NULL;
  mrc_obj_get_param_string((struct mrc_obj*)out, "data_dir" ,&outdir);

  mrc_io_set_name(out_c->io, "particles_out");
  mrc_io_set_param_string(out_c->io, "basename", "particles");
  mrc_io_set_param_string(out_c->io, "outdir",outdir);
  mrc_io_setup(out_c->io);
  mrc_io_view(out_c->io);
  
  out_c->next=out_c->first;
}

static void
psc_output_particles_custom_binary_destroy(struct psc_output_particles *out)
{
  struct psc_output_particles_custom_binary *out_c = to_psc_output_particles_custom_binary(out);
  mrc_io_destroy(out_c->io);
}

static void
psc_output_particles_custom_binary_set_from_options(struct psc_output_particles *out)
{
  struct psc_output_particles_custom_binary *out_c = to_psc_output_particles_custom_binary(out);
  mrc_io_set_from_options(out_c->io);
}

#define VAR(x) (void *)offsetof(struct psc_output_particles_custom_binary, x)

static struct param psc_output_particles_custom_binary_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },//This default value will NOT be used. The defaul value of the parameter outdir of mrc_io will be used instead.
  { "filter_function"   , VAR(filter_func) , PARAM_PTR(NULL) },
  { "particle_first"       , VAR(first)         , PARAM_INT(0)            },
  { "particle_step"        , VAR(step)          , PARAM_INT(10)           },
  { "write_charge"        , VAR(write_q)          , PARAM_BOOL(true)           },
  { "write_mass"        , VAR(write_m)          , PARAM_BOOL(true)           },
  { "write_x"        , VAR(write_x)          , PARAM_BOOL(true)           },
  { "write_y"        , VAR(write_y)          , PARAM_BOOL(true)           },
  { "write_z"        , VAR(write_z)          , PARAM_BOOL(true)           },
  { "write_px"        , VAR(write_px)          , PARAM_BOOL(true)           },
  { "write_py"        , VAR(write_py)          , PARAM_BOOL(true)           },
  { "write_pz"        , VAR(write_pz)          , PARAM_BOOL(true)           },
  {},
};
#undef VAR



// ======================================================================
// psc_output_particles: subclass "c"

struct psc_output_particles_ops psc_output_particles_custom_binary_ops = {
  .name                  = "custom_binary",
  .size                  = sizeof(struct psc_output_particles_custom_binary),
  .param_descr           = psc_output_particles_custom_binary_descr,
  .create                = psc_output_particles_custom_binary_create,
  .setup                 = psc_output_particles_custom_binary_setup,
  .set_from_options      = psc_output_particles_custom_binary_set_from_options,
  .destroy               = psc_output_particles_custom_binary_destroy,
//  .write                 = psc_output_particles_c_write,
//  .read                  = psc_output_particles_c_read,
//  .view                  = psc_output_particles_c_view,
  .run                   = psc_output_particles_custom_binary_run,
};
