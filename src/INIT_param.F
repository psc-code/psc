
#include "config.h"

c     ----------------------------------------------------------------------
      subroutine INIT_basic
c     ----------------------------------------------------------------------

      use VLA_variables

      implicit none
      include 'mpif.h'

c PARAMETERS THAT SHOULD NOT BE CHANGED UNDER NORMAL CONDITIONS

c pi: pi

      pi=3.1415927

      call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,info)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mpe,info)

      end subroutine INIT_basic

c PARAMETER INITIALIZATION 

c     ----------------------------------------------------------------------
      subroutine INIT_param
c     ----------------------------------------------------------------------

      use VLA_variables
      use PSC_globals
      use PSC_domain
      use PSC_patch

      implicit none

c REQUIRED SETTINGS
c Please do not change


      if (i1x==i1n) then
         xnpe=1
         boundary_part_x = 1
         boundary_fld_x1 = 1
         boundary_fld_x2 = 1
      else
         if (i1x-i1n+1.le.(rd1+1)*xnpe) then
            i1x=i1n+(rd1+1)*xnpe-1
         endif
         if (boundary_fld_x1==2.or.boundary_fld_x2==2) then
            if (i1x-i1n+1.le.2*(size+1)) then 
               i1x=i1n+2*(size+1)-1
            endif
            if (i1x-i1n+1.le.(size+1)*xnpe) then
               i1x=i1n+(size+1)*xnpe-1
            endif
         endif
         if (boundary_fld_x1==3.or.boundary_fld_x2==3) then
            if (i1x-i1n+1.le.2*(size+1)) then 
               i1x=i1n+2*(size+1)-1
            endif
            if (i1x-i1n+1.le.(size+1)*xnpe) then
               i1x=i1n+(size+1)*xnpe-1
            endif
         endif
         if (i1x-i1n+1.gt.i1tot) then
            i1tot=i1x-i1n+1
         endif
      endif
      if (i2x==i2n) then
         ynpe=1
         boundary_part_y = 1
         boundary_fld_y1 = 1
         boundary_fld_y2 = 1
      else
         if (i2x-i2n+1.le.(rd2+1)*ynpe) then
            i2x=i2n+(rd2+1)*ynpe-1
         endif
         if (boundary_fld_y1==2.or.boundary_fld_y2==2) then
            if (i2x-i2n+1.le.2*(size+1)) then  
               i2x=i2n+2*(size+1)-1
            endif
            if (i2x-i2n+1.le.(size+1)*ynpe) then
               i2x=i2n+(size+1)*ynpe-1
            endif
         endif
         if (boundary_fld_y1==3.or.boundary_fld_y2==3) then
            if (i2x-i2n+1.le.2*(size+1)) then  
               i2x=i2n+2*(size+1)-1
            endif
            if (i2x-i2n+1.le.(size+1)*ynpe) then
               i2x=i2n+(size+1)*ynpe-1
            endif
         endif
         if (i2x-i2n+1.gt.i2tot) then
            i2tot=i2x-i2n+1
         endif
      endif
      if (i3x==i3n) then
         znpe=1
         boundary_part_z = 1
         boundary_fld_z1 = 1
         boundary_fld_z2 = 1
      else
         if (i3x-i3n+1.le.(rd3+1)*znpe) then
            i3x=i3n+(rd3+1)*znpe-1
         endif
         if (boundary_fld_z1==2.or.boundary_fld_z2==2) then
            if (i3x-i3n+1.le.2*(size+1)) then
               i3x=i3n+2*(size+1)-1
            endif
            if (i3x-i3n+1.le.(size+1)*znpe) then
               i3x=i3n+(size+1)*znpe-1
            endif
         endif
         if (boundary_fld_z1==3.or.boundary_fld_z2==3) then
            if (i3x-i3n+1.le.2*(size+1)) then
               i3x=i3n+2*(size+1)-1
            endif
            if (i3x-i3n+1.le.(size+1)*znpe) then
               i3x=i3n+(size+1)*znpe-1
            endif
         endif
         if (i3x-i3n+1.gt.i3tot) then
            i3tot=i3x-i3n+1
         endif
      endif


c wl: laser frequency
c phi0: electric potential normalization
c a0: vector potential normalization
c e0: electric field normalization
c b0: magnetic field normalization
c j0: current density normalization
c rho0: charge density normalization
c ld: laser wave number
c vos: oscillation velocity
c vt: thermal velocity
c wp: plasma frequency


      wl=2.0*pi*cc/lw
      ld=cc/wl
      if (e0 == 0.) then
         e0=dsqrt(2.0*i0/eps0/cc)/lw/1.0e6
      endif
      b0=e0/cc
      j0=eps0*wl*e0
      rho0=eps0*wl*b0
      phi0=ld*e0
      a0=e0/wl
      vos=qq*e0/(mm*wl)
      vt=dsqrt(tt/mm)
      wp=dsqrt(qq**2*n0/eps0/mm)
      alpha=wp/wl
      beta=vt/cc
      eta=vos/cc


c dx: spatial resolution in x
c dy: spatial resolution in y
c dz: spatial resolution in z
c dt: time step, Courant criterion
c nnp: number of timesteps for a full laser cycle
c np:  number of timesteps for time-averaging

      ! KG FIXME, is i[123]tot(not -1) right in all cases?
      dx=lengthx/(i1tot)/ld
      dy=lengthy/(i2tot)/ld
      dz=lengthz/(i3tot)/ld

      dt=0.75*sqrt(1.0/(1.0/dx**2+1.0/dy**2+1.0/dz**2))

c     adjust to match laser cycles
      nnp=int(2.0*pi/dt)+1
      dt=2.0*pi/nnp
      np=2*nnp
      if (nmax == 0) then
         nmax=30*nnp
      endif

      end subroutine INIT_param

c     ----------------------------------------------------------------------
      subroutine INIT_param_fortran
c     ----------------------------------------------------------------------
c     setup processor map

      use VLA_variables

      implicit none

      if (use_pml) call PML_coeff_init_param

! position condition for time dependent pml

      pos_x1=0.0
      pos_x2=0.0
      pos_y1=0.0
      pos_y2=0.0
      pos_z1=0.0
      pos_z2=0.0


c I/O CONTROL PARAMETERS

c data_out: path of data output directory
c data_chk: path of checkpointing directory

      data_out='.'
      data_chk='.'
      data_processed='.'
      
c OTHER SETTINGS THAT SHOULD NOT BE CHANGED
c UNDER NORMAL CONDITIONS

c nstart: time step counter for starting the code
c shift_c: time step counter for shifting the frame in the code
c shift_z: shifts the density function in the code
c cpuinou: WALL CLOCK time required for I/O
c cpucomp: WALL CLOCK time required for computation
c cpus: WALL CLOCK time at end time


      nstart=0
      shift_c=0
      shift_l=0
      shift_z=0.0

      end subroutine INIT_param_fortran

c     ----------------------------------------------------------------------
      subroutine INIT_grid_map
c     ----------------------------------------------------------------------

      use PSC_domain

      implicit none

      integer :: i1, i2, i3, pec

      if (allocated(seg_i1)) then
         return
      endif

      allocate(seg_i1(0:xnpe*ynpe*znpe-1))
      allocate(seg_i2(0:xnpe*ynpe*znpe-1))
      allocate(seg_i3(0:xnpe*ynpe*znpe-1))
      allocate(seg_inv(1:xnpe,1:ynpe,1:znpe))

      do i3=1,znpe
         do i2=1,ynpe
            do i1=1,xnpe
               pec=(i3-1)*xnpe*ynpe+(i2-1)*xnpe+(i1-1)
               seg_i1(pec)=i1
               seg_i2(pec)=i2
               seg_i3(pec)=i3
               seg_inv(i1,i2,i3)=pec
            enddo
         enddo
      enddo

      end subroutine INIT_grid_map


