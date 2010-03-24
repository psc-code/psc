# -*- Makefile -*- GNU-Make only (what else ?)
#
# Author:   Hartmut Ruhl
# File:     Makefile for PSC code
# Date:     10/03/2005
#
#################################################


# COMPILERS AND FLAGS TO USE


############################################################################################################
# SETTINGS FOR THE CRAY T3E

#FFLAGS= -R ab -g
#FFLAGS= -O3,unroll2,inline4,pipeline3
#LIB= -lmfastv 
#FC= f90

#FFLAGS= -fast -check bounds -check format -warn unused -I/usr/local/include/
#FC= mpif90

#FFLAGS= -fast -warn unused
#FC= mpif90

#FFLAGS= -O4
#FC= mpxlf
############################################################################################################


############################################################################################################
# SETTINGS FOR THE NERSC SEABORG

#FFLAGS= -O2 -qstrict -qarch = pwr3 -bmaxdata:0x70000000 -bmaxstack:0x10000000 -qfixed -I../Includes
#FC= mpxlf90
############################################################################################################


############################################################################################################
# SETTINGS FOR THE EARTH SIMULATOR

#FFLAGS= -Wf"-pvctl loopcnt=335544320" -C hopt -ftrace
#FC=esmpif90
#AR=esar
############################################################################################################


############################################################################################################
# SETTINGS FOR PICBOX

#FFLAGS= -O3 -parallel  -I/usr/local/include/
#FFLAGS=-I/usr/local/include/ -O3
#LIB= -lmpich -lfmpich -lpthread
#FC=ifort
############################################################################################################


############################################################################################################
# SETTINGS FOR LINUX OPTERONS

#FFLAGS= -Mbounds
#FFLAGS= -fast
#FC=mpif90
FFLAGS= -O3

#FC=/opt/mpi-gm-gfortran/bin/mpif90
# CC=g++
# CFLAGS= -O3
# LIB=-lstdc++
############################################################################################################


############################################################################################################
# SETTINGS FOR MAC OS X

FFLAGS= -Rb -Rc -Rs
FFLAGS= -O3
FC=mpif90
############################################################################################################



# Other Vars


BASE=$(shell basename `pwd`)


# Module Objects

IOBJS= VLA_variables.o PIC_variables.o VLI.o \
INIT_ppulse_x1.o INIT_spulse_x1.o \
INIT_ppulse_x2.o INIT_spulse_x2.o \
INIT_ppulse_y1.o INIT_spulse_y1.o \
INIT_ppulse_y2.o INIT_spulse_y2.o \
INIT_ppulse_z1.o INIT_spulse_z1.o \
INIT_ppulse_z2.o INIT_spulse_z2.o \
PIC_bpulse.o PIC_ionize.o \
INIT_field.o INIT_param.o INIT_den.o INIT_idistr.o \
PIC_pml_msa.o PIC_pml_msb.o \
PIC_msa.o PIC_msb.o \
PIC_fex.o PIC_fey.o PIC_fez.o \
PIC_fax.o PIC_fay.o PIC_faz.o \
SERV_labelgen.o SERV_systime.o \
SERV_write.o \
OUT_param.o OUT_field.o OUT_part.o OUT_poyc.o \
OUT_count.o OUT_energy.o PIC_pex.o PIC_pey.o PIC_pez.o \
PIC_move_part.o PIC_move_part_xyz.o PIC_move_part_xy.o \
PIC_move_part_xz.o PIC_move_part_yz.o PIC_move_part_x.o \
PIC_move_part_y.o PIC_move_part_z.o PIC_bin_coll.o PIC_sort.o \
MCC_variables.o MCC_init.o MCC_impact.o


MOBJS= VLA_variables.o PIC_variables.o VLA.o \
INIT_ppulse_x1.o INIT_spulse_x1.o \
INIT_ppulse_x2.o INIT_spulse_x2.o \
INIT_ppulse_y1.o INIT_spulse_y1.o \
INIT_ppulse_y2.o INIT_spulse_y2.o \
INIT_ppulse_z1.o INIT_spulse_z1.o \
INIT_ppulse_z2.o INIT_spulse_z2.o \
PIC_bpulse.o PIC_ionize.o \
INIT_field.o INIT_param.o INIT_den.o INIT_idistr.o \
PIC_pml_msa.o PIC_pml_msb.o \
PIC_msa.o PIC_msb.o \
PIC_fex.o PIC_fey.o PIC_fez.o \
PIC_fax.o PIC_fay.o PIC_faz.o \
SERV_labelgen.o SERV_systime.o \
SERV_write.o SERV_read.o \
OUT_param.o OUT_field.o OUT_part.o OUT_poyc.o \
OUT_count.o OUT_energy.o PIC_pex.o PIC_pey.o PIC_pez.o \
PIC_move_part.o PIC_move_part_xyz.o PIC_move_part_xy.o \
PIC_move_part_xz.o PIC_move_part_yz.o PIC_move_part_x.o \
PIC_move_part_y.o PIC_move_part_z.o PIC_bin_coll.o PIC_sort.o \
MCC_variables.o MCC_init.o MCC_impact.o

SOBJS= VLA_variables.o PIC_variables.o INIT_param.o \
SELECT.o SERV_labelgen.o \
SELECT_pfield_evol.o SELECT_tfield_evol.o \
SELECT_cl_evol.o SELECT_count_evol.o \
SELECT_electron_evol.o SELECT_ion_evol.o SELECT_atom_evol.o

# Main Targets

ALL: VLI.x VLA.x SELECT.x BACKUP

# uncomment for cpp
VLI.x: libPICI.a
	$(FC) -o $@ $(FFLAGS) VLI.o libPICI.a $(LIB)
VLA.x: libPICM.a
	$(FC) -o $@ $(FFLAGS) VLA.o libPICM.a $(LIB)
# VLI.x: libPICI.a
# 	$(FC) -o $@ $(FFLAGS) VLI.o libPICI.a
# VLA.x: libPICM.a
# 	$(FC) -o $@ $(FFLAGS) VLA.o libPICM.a
SELECT.x: libPICS.a 
	$(FC) -o $@ $(FFLAGS) SELECT.o libPICS.a

libPICI.a: $(IOBJS)
	$(AR) -ruvs libPICI.a $(IOBJS)
libPICM.a: $(MOBJS)
	$(AR) -ruvs libPICM.a $(MOBJS)
libPICS.a: $(SOBJS)
	$(AR) -ruvs libPICS.a $(SOBJS)


# Special Dependencies


INIT_spulse_x1.o: INIT_spulse_x1.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_spulse_x2.o: INIT_spulse_x2.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_spulse_y1.o: INIT_spulse_y1.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_spulse_y2.o: INIT_spulse_y2.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_spulse_z1.o: INIT_spulse_z1.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_spulse_z2.o: INIT_spulse_z2.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_ppulse_x1.o: INIT_ppulse_x1.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_ppulse_x2.o: INIT_ppulse_x2.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_ppulse_y1.o: INIT_ppulse_y1.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_ppulse_y2.o: INIT_ppulse_y2.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_ppulse_z1.o: INIT_ppulse_z1.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_ppulse_z2.o: INIT_ppulse_z2.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_bpulse.o: PIC_bpulse.f PIC_variables.o VLA_variables.o SERV_systime.o
	$(FC) $(FFLAGS) -c $<

PIC_ionize.o: PIC_ionize.f PIC_variables.o VLA_variables.o SERV_systime.o
	$(FC) $(FFLAGS) -c $<

INIT_field.o: INIT_field.f PIC_fax.o PIC_fay.o PIC_faz.o \
              PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_idistr.o: INIT_idistr.f SERV_labelgen.o PIC_variables.o \
	       VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_param.o: INIT_param.f PIC_variables.o \
	      VLA_variables.o
	$(FC) $(FFLAGS) -c $<

INIT_den.o: INIT_den.f
	$(FC) $(FFLAGS) -c $<

OUT_field.o: OUT_field.f SERV_systime.o SERV_labelgen.o \
	     PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

OUT_param.o: OUT_param.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

OUT_part.o: OUT_part.f SERV_systime.o SERV_labelgen.o \
	    PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

OUT_energy.o: OUT_energy.f SERV_systime.o SERV_labelgen.o \
	    PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

OUT_poyc.o: OUT_poyc.f SERV_systime.o SERV_labelgen.o \
	    PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_bin_coll.o: PIC_bin_coll.f SERV_systime.o SERV_labelgen.o \
	        PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

MCC_init.o: MCC_init.f PIC_variables.o VLA_variables.o MCC_variables.o
	$(FC) $(FFLAGS) -c $<
MCC_impact.o: MCC_impact.f PIC_variables.o VLA_variables.o MCC_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_pex.o: PIC_pex.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_fex.o: PIC_fex.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_fax.o: PIC_fax.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_pey.o: PIC_pey.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_fey.o: PIC_fey.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_fay.o: PIC_fay.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_pez.o: PIC_pez.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_fez.o: PIC_fez.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_faz.o: PIC_faz.f PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_move_part.o: PIC_move_part.f PIC_move_part_xyz.o PIC_move_part_xy.o PIC_move_part_xz.o \
		 PIC_move_part_yz.o PIC_move_part_x.o PIC_move_part_y.o PIC_move_part_z.o \
		 PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_move_part_xyz.o: PIC_move_part_xyz.f SERV_systime.o PIC_fax.o PIC_fay.o PIC_faz.o \
		PIC_pex.o PIC_pey.o PIC_pez.o \
		PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_move_part_xy.o: PIC_move_part_xy.f SERV_systime.o PIC_fax.o PIC_fay.o \
		 PIC_pex.o PIC_pey.o \
		 PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_move_part_xz.o: PIC_move_part_xz.f SERV_systime.o PIC_fax.o PIC_faz.o \
                 PIC_pex.o PIC_pez.o \
		 PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_move_part_yz.o: PIC_move_part_yz.f SERV_systime.o PIC_fay.o PIC_faz.o \
		PIC_pey.o PIC_pez.o \
		PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_move_part_x.o: PIC_move_part_x.f SERV_systime.o PIC_fax.o PIC_pex.o \
		 PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_move_part_y.o: PIC_move_part_y.f SERV_systime.o PIC_fay.o PIC_pey.o \
		 PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<
PIC_move_part_z.o: PIC_move_part_z.f SERV_systime.o PIC_faz.o PIC_pez.o \
		 PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_pml_msa.o: PIC_pml_msa.f SERV_systime.o  PIC_fex.o PIC_fey.o PIC_fez.o \
           PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_msa.o: PIC_msa.f SERV_systime.o  PIC_fex.o PIC_fey.o PIC_fez.o \
           PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_psa.o: PIC_psa.f SERV_systime.o  PIC_fex.o PIC_fey.o PIC_fez.o \
           PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_pml_msb.o: PIC_pml_msb.f SERV_systime.o  PIC_fex.o PIC_fey.o PIC_fez.o \
           PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_msb.o: PIC_msb.f SERV_systime.o  PIC_fex.o PIC_fey.o PIC_fez.o \
           PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

PIC_sort.o: PIC_sort.f SERV_systime.o PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SELECT.o: SELECT.f INIT_param.o SELECT_pfield_evol.o SELECT_tfield_evol.o \
	  SELECT_cl_evol.o SELECT_electron_evol.o SELECT_ion_evol.o \
          SELECT_atom_evol.o PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SELECT_atom_evol.o: SELECT_atom_evol.f SERV_labelgen.o \
		    PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SELECT_cl_evol.o: SELECT_cl_evol.f SERV_labelgen.o \
                  PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SELECT_count_evol.o: SELECT_count_evol.f SERV_labelgen.o \
                  PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SELECT_electron_evol.o: SELECT_electron_evol.f SERV_labelgen.o \
			PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SELECT_ion_evol.o: SELECT_ion_evol.f SERV_labelgen.o \
		   PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SELECT_pfield_evol.o: SELECT_pfield_evol.f SERV_labelgen.o \
		      PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SELECT_tfield_evol.o: SELECT_tfield_evol.f SERV_labelgen.o \
		      PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SERV_read.o: SERV_read.f SERV_labelgen.o PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

SERV_write.o: SERV_write.f SERV_labelgen.o PIC_variables.o VLA_variables.o
	$(FC) $(FFLAGS) -c $<

VLI.o: VLI.f INIT_param.o SERV_systime.o INIT_idistr.o INIT_field.o \
       OUT_param.o PIC_sort.o PIC_bin_coll.o \
       PIC_pml_msa.o PIC_pml_msb.o OUT_field.o OUT_part.o \
       PIC_msa.o PIC_msb.o \
       OUT_poyc.o PIC_move_part.o SERV_write.o \
       PIC_variables.o VLA_variables.o SERV_consist.f \
       SERV_cput.f OUT_cput.f
	$(FC) $(FFLAGS) -c $<

VLA.o: VLA.f INIT_param.o SERV_systime.o INIT_idistr.o INIT_field.o \
       OUT_param.o PIC_sort.o PIC_bin_coll.o \
       PIC_pml_msa.o PIC_pml_msb.o OUT_field.o OUT_part.o \
       PIC_msa.o PIC_msb.o \
       OUT_poyc.o PIC_move_part.o SERV_write.o SERV_read.o \
       PIC_variables.o VLA_variables.o SERV_consist.f \
       SERV_cput.f OUT_cput.f
	$(FC) $(FFLAGS) -c $<


# Default Rules


%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# %.o: %.cpp
# 	$(CC) $(CFLAGS) -c $<

clean: cleano cleanx cleana cleanmod cleaneo cleangz

cleano:
	rm -f PP TT EE II AA RMFILES *FILES core
	rm -f okfile endfile *out *~ *.gz *.eps *.data
	rm -f PI1* PI2* PI3* PI4* PI5* PI6* PI7* PI8* PI9*

cleanx: 
	rm -f *.x

cleana:
	rm -f *.a

cleanmod:
	rm -f *.mod

cleaneo:
	rm -f *.e* *.o*

cleangz:
	rm -f *.gz


BACKUP: Makefile *.f *.cpp *.pro *exec PROCESSOR_* vl*exec IDL* idl_* LN LS MKDIR RM DU KILL Machines mkanim
	( tar cvf `date +%Y.%m.%d`.tar ./Makefile ./*.f *.cpp \
         ./*.pro ./*exec ./PROCESSOR_* ./vl*exec ./IDL* ./idl_* LN LS MKDIR RM DU KILL Machines mkanim; \
         mv ./*tar ./BACKUP )
