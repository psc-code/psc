
#include "config.h"

      module VLA_variables
      implicit none

#ifdef USE_PML
      real(kind=8) dx2A,dy2A,dz2A,hx2A,hy2A,hz2A  ! ab  
      real(kind=8) dx2B,dy2B,dz2B,hx2B,hy2B,hz2B  ! ab
      real(kind=8) poyx1i,poyx1o,poyx2i,poyx2o                 ! ab
      real(kind=8) poyy1i,poyy1o,poyy2i,poyy2o                 ! ab
      real(kind=8) poyz1i,poyz1o,poyz2i,poyz2o                 ! ab
      real(kind=8) poyni, poyno                 ! ab
#endif

c seg_i1: x-label of distributed data segment
c seg_i2: y-label of distributed data segment
c seg_i3: z-label of distributed data segment
c seg_inv: returns node number after x-, y-label, and z-label input

      integer,allocatable, dimension(:) :: seg_i1,seg_i2,seg_i3
      integer,allocatable, dimension(:,:,:) :: seg_inv

c cpumess: WALL CLOCK time for cumulative communication
c cpucomp: WALL CLOCK time for cumulative computation
c cpuinou: WALL CLOCK time for cumulative input/output
c cpus: WALL CLOCK time at program start
c cpuf: WALL CLOCK time at program stop
c cpue: MAXIMUM WALL CLOCK time of all PEs at program end
c cpum: WALL CLOCK time limit
c cpu: stores WALL CLOCK times on individual PEs
c cpu_ary: stores WALL CLOCK times of all PEs

      real(kind=8) cpumess,cpucomp,cpuinou
      real(kind=8) cpua,cpub
      real(kind=8) cpuc,cpud
      real(kind=8) cpug,cpuh
      real(kind=8) cpus,cpuf,cpue,cpum
      real(kind=8),allocatable,dimension(:) :: cpu
      real(kind=8),allocatable,dimension(:,:) :: cpu_ary

c nprf: time step counter for data output
c nprc: time step counter for collision data output
c tmnvf: time step counter for starting time averaging of fields
c tmxvf: time step counter for ending time averaging of fields
c tmnvp: time step counter for starting time averaging of poynting flux
c tmxvp: time step counter for ending time averaging of poynting flux
c tmnvc: time step counter for starting time averaging of particle counting
c tmxvc: time step counter for ending time averaging of particle counting
c nmax: maximum number of time steps 
c nstart: time step counter for restarting the code
c shift_c: time step counter for moving the frame in the code
c shift_l: step counter for moving the frame in the code
c nnp: number of timesteps for one full cycle
c np:  number of timesteps for time-averaging

      integer nprfo,nprf,nprco,nprc
      integer tmnvfo,tmxvfo,tmnvf,tmxvf
      integer tmnvpo,tmxvpo,tmnvp,tmxvp
      integer tmnvco,tmxvco,tmnvc,tmxvc
      integer nmax,nstart,shift_c,shift_l,np,nnp

      real(kind=8) fluxi,fluxo
      real(kind=8) ex2A,ey2A,ez2A,bx2A,by2A,bz2A
      real(kind=8) ex2B,ey2B,ez2B,bx2B,by2B,bz2B
      real(kind=8) pox,poy,poz,je
      real(kind=8) p2A,p2B

c i0: irradiance in Wum**2/m**2
c n0: peak atom density in m**3
c lw: laser wavelength
c wl: laser frequency
c wp: plasma frequency
c vt: thermal velocity
c ld: plasma wavelength
c qq: reference charge in As
c mm: reference mass in kg
c tt: atom temperature in J
c cc: light velocity in m/s
c eps0: eps0 in As/Vm
c pi: the number pi
c dt: time resolution in dimensionless units
c dx: gridspacing in x-direction in dimensionless units
c dy: gridspacing in y-direction in dimensionless units
c dz: gridspacing in z-direction in dimensionless units
c shift_z: density function shift in dimensionless units
c xnpe: number of PEs in x-direction
c ynpe: number of PEs in y-direction
c znpe: number of PEs in z-direction
c i1mn: local size of the grid in x on a slave PE
c i1mx: local size of the grid in x on a slave PE
c i2mn: local size of the grid in y on a slave PE
c i2mx: local size of the grid in y on a slave PE
c i3mn: local size of the grid in z on a slave PE
c i3mx: local size of the grid in z on a slave PE
c npe: total number of PEs in the partition
c mpe: number of the PE the code is running on
c dnprf: time increment for output of fields
c dnprc: time increment for output of collision data
c lengthx: x extension of simulation box in m
c lengthy: y extension of simulation box in m
c lengthz: z extension of simulation box in m
c info: error message from mpi routines
c pec: PE counter
c r1n: min output range along x in grid points
c r1x: max output range along x in grid points
c r2n: min output range along y in grid points
c r2x: max output range along y in grid points
c r3n: min output range along x in grid points
c r3x: max output range along x in grid points

! deltax,y,z: width of pml in dimensionless units

      real(kind=8) i0,n0,lw,wl
      real(kind=8) qq,mm,tt,cc,eps0,pi
      real(kind=8) ld,wp,vos,vt
      real(kind=8) dt,dx,dy,dz,shift_z
      real(kind=8) alpha,beta,eta
      real(kind=8) e0,b0,j0,rho0, phi0, a0
      real(kind=8) lengthx,lengthy,lengthz
      real(kind=8) deltax, deltay, deltaz            ! added by ab

      integer info,pec
      integer i1,i2,i3
      integer n,i,j,k,l
      integer i1tot,i2tot,i3tot
      integer i1n,i1x,i2n,i2x,i3n,i3x
      integer r1n,r1x,r2n,r2x,r3n,r3x

      integer xnpe,ynpe,znpe
      integer i1mn,i1mx,i2mn,i2mx,i3mn,i3mx
      integer npe,mpe
      integer rd1,rd2,rd3,dnprf,dnprc

c data_out: data out directory
c data_chk: data checkpointing directory

      character*(200) data_out
      character*(200) data_chk
      character*(20) char
      character*(5) pe

c boundary_part_x: particle boundary conditions in x
c boundary_part_y: particle boundary conditions in y
c boundary_part_z: particle boundary conditions in z

      integer boundary_part_x
      integer boundary_part_y
      integer boundary_part_z

! boundary_fld_x1, boundary_fld_x2: boundary condition for bnd in x-front,back
! boundary_fld_y1, boundary_fld_y2: boundary condition for bnd in y-left,right
! boundary_fld_z1, boundary_fld_z2: boundary condition for bnd in z-up,down
!
! condition = 0 : radiating boundaries
! condition = 1 : periodic boundaries
! condition = 2 : upmls
! condition = 3 : radiating boundaries with upmls at a given time

      integer boundary_fld_x1
      integer boundary_fld_x2
      integer boundary_fld_y1
      integer boundary_fld_y2
      integer boundary_fld_z1
      integer boundary_fld_z2

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
      subroutine PIC_set_variables(
     *     i1mn_, i2mn_, i3mn_,
     *     i1mx_, i2mx_, i3mx_,
     *     rd1_, rd2_, rd3_,
     *     i1n_, i2n_, i3n_,
     *     i1x_, i2x_, i3x_,
     *     cori_, alpha_, eta_,
     *     dt_, dx_, dy_, dz_,
     *     wl_, wp_, n_)
c     ======================================================================

      use PIC_variables
      use VLA_variables

      implicit none

      integer :: i1mn_, i2mn_, i3mn_, i1mx_, i2mx_, i3mx_
      integer :: rd1_, rd2_, rd3_
      integer :: i1n_, i2n_, i3n_, i1x_, i2x_, i3x_
      integer :: n_
      real(kind=8) :: cori_, alpha_, eta_
      real(kind=8) :: dt_, dx_, dy_, dz_
      real(kind=8) :: wl_, wp_

      i1mn = i1mn_
      i2mn = i2mn_
      i3mn = i3mn_
      i1mx = i1mx_
      i2mx = i2mx_
      i3mx = i3mx_
      rd1 = rd1_
      rd2 = rd2_
      rd3 = rd3_
      i1n = i1n_
      i2n = i2n_
      i3n = i3n_
      i1x = i1x_
      i2x = i2x_
      i3x = i3x_
      
      cori = cori_
      alpha = alpha_
      eta = eta_
      dt = dt_
      dx = dx_
      dy = dy_
      dz = dz_

      wl = wl_
      wp = wp_
      n = n_

      end subroutine PIC_set_variables

c     ======================================================================
      module fields
c     ======================================================================
      implicit none

c ex-bz: em field arrays 
c ne-nn: charge densities
c jxi-jzi: current densities

      real(kind=8),allocatable,dimension(:,:,:) :: ex,ey,ez
      real(kind=8),allocatable,dimension(:,:,:) :: bx,by,bz
      real(kind=8),allocatable,dimension(:,:,:) :: jxi,jyi,jzi
      real(kind=8),allocatable,dimension(:,:,:) :: ne,ni,nn

      end module fields

