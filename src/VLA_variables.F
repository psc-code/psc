
#include "config.h"

      module VLA_variables
      implicit none

      logical use_pml

      real(kind=8) dx2A,dy2A,dz2A,bx2A,by2A,bz2A  ! ab  
      real(kind=8) dx2B,dy2B,dz2B,bx2B,by2B,bz2B  ! ab
      real(kind=8) poyx1i,poyx1o,poyx2i,poyx2o                 ! ab
      real(kind=8) poyy1i,poyy1o,poyy2i,poyy2o                 ! ab
      real(kind=8) poyz1i,poyz1o,poyz2i,poyz2o                 ! ab
      real(kind=8) poyni, poyno                 ! ab

c nmax: maximum number of time steps 
c nstart: time step counter for restarting the code
c shift_c: time step counter for moving the frame in the code
c shift_l: step counter for moving the frame in the code
c nnp: number of timesteps for one full cycle
c np:  number of timesteps for time-averaging


      integer nmax,nstart,shift_c,shift_l,np,nnp

      real(kind=8) fluxi,fluxo
      real(kind=8) ex2A,ey2A,ez2A,hx2A,hy2A,hz2A
      real(kind=8) ex2B,ey2B,ez2B,hx2B,hy2B,hz2B
      real(kind=8) pox,poy,poz,je
      real(kind=8) p2A,p2B

c pi: the number pi
c shift_z: density function shift in dimensionless units
c npe: total number of PEs in the partition
c mpe: number of the PE the code is running on
c lengthx: x extension of simulation box in m
c lengthy: y extension of simulation box in m
c lengthz: z extension of simulation box in m
c info: error message from mpi routines
! deltax,y,z: width of pml in dimensionless units

      real(kind=8) pi
      real(kind=8) shift_z
      real(kind=8) lengthx,lengthy,lengthz
      real(kind=8) deltax, deltay, deltaz            ! added by ab

      integer info
      integer n

      integer npe,mpe

c data_out: data out directory
c data_chk: data checkpointing directory

      character*(200) data_out
      character*(200) data_chk
      character*(200) data_processed
      character*(20) char
      character*(5) pe

! checking condition for time dependent pml

      real(kind=8) pos_x1
      real(kind=8) pos_x2
      real(kind=8) pos_y1
      real(kind=8) pos_y2
      real(kind=8) pos_z1
      real(kind=8) pos_z2

! thick: thickness of pml in gridpoints
! size: thickness of pml and buffer region in gridpoints

      integer thick, size

! pmlcheck hacky globals...

      integer :: pml_xmin,pml_xmax     ! ab
      integer :: pml_ymin,pml_ymax     ! ab
      integer :: pml_zmin,pml_zmax     ! ab

      end module VLA_variables

c     ======================================================================
      subroutine PIC_set_variables(wl_, wp_, n_)
c     ======================================================================

      use VLA_variables, only: n
      use PSC_params
      implicit none

      integer :: n_
      real(kind=8) :: wl_, wp_

      wl = wl_
      wp = wp_
      n = n_

      end subroutine PIC_set_variables

