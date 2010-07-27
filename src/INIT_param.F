
#include "config.h"

c     ----------------------------------------------------------------------
      subroutine INIT_basic
c     ----------------------------------------------------------------------

      use VLA_variables
      use PIC_variables

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
      use PIC_variables

      implicit none


c THESE PARAMETRES ARE USED FOR DETERMINING GRID SIZE,
c PARTICLE NUMBER, DATA PARTITIONING, AND PHYSICAL
c PARAMETERS!

#ifdef USE_C_INTERFACE
      call C_INIT_param
#else
      call INIT_param_domain
      call INIT_param_psc
#endif

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


c QUASI-PARTICLE MASS FACTOR

#ifndef USE_C_INTERFACE
      if (nicell.gt.0) then
         cori=1.0/nicell
      else
         cori=1.0e20
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
#endif

      end subroutine INIT_param

c     ----------------------------------------------------------------------
      subroutine INIT_param_fortran
c     ----------------------------------------------------------------------
c     setup processor map

      use VLA_variables
      use PIC_variables

      implicit none

! PML PARAMETERS - added by ab

! deltax: thickness of pml in x direction
! deltay: thickness of pml in y direction
! deltaz: thickness of pml in z direction

      deltax = thick*dx
      deltay = thick*dy
      deltaz = thick*dz

! attenuation factors
      
      kappax_max = 1.34
      sigmax_max = 13.8*eps0/deltax

      kappay_max = 1.34
      sigmay_max = 13.8*eps0/deltay

      kappaz_max = 1.34
      sigmaz_max = 13.8*eps0/deltaz

! position condition for time dependent pml

      pos_x1=0.0
      pos_x2=0.0
      pos_y1=0.0
      pos_y2=0.0
      pos_z1=0.0
      pos_z2=0.0


c I/O CONTROL PARAMETERS

c nprf: first time step for output of fields
c dnprf: timestep increment for output of fields
c r1n: min output range along x
c r1x: max output range along x
c r2n: min output range along y
c r2x: max output range along y
c r3n: min output range along z
c r3x: max output range along z

      nprf=0
      dnprf=10                                      ! changed by ab
      r1n=i1n
      r1x=i1x
      r2n=i2n
      r2x=i2x
      r3n=i3n
      r3x=i3x

c nprc: first time step for output of collision parameters
c dnprc: timestep increment for output of collision parameters

      nprc=0
      dnprc=100

c nprparti: first time step for output of ions
c dnprparti: timestep increment for output of ions
c nistep: index increment for ion output
c plin: smallest particle number for output
c plix: largest particle number for output

      nprparti=0
      dnprparti=10
      nistep=1
      plin=0
      plix=30000000

c tmnvf: starting time averaging of fields
c tmxvf: ending time averaging of fields
c tmnvp: starting time averaging of poynting flux
c tmxvp: ending time averaging of poynting flux
c tmnvc: starting time averaging of counting
c tmxvc: ending time averaging of counting

      tmnvf=0*nnp+1
      tmxvf=0*nnp+np
      tmnvp=0*nnp+1         
      tmxvp=0*nnp+np
      tmnvc=0*nnp+1         
      tmxvc=0*nnp+np

c data_out: path of data output directory
c data_chk: path of checkpointing directory

      data_out='.'
      data_chk='.'

c OTHER SETTINGS THAT SHOULD NOT BE CHANGED
c UNDER NORMAL CONDITIONS

c nstart: time step counter for starting the code
c shift_c: time step counter for shifting the frame in the code
c shift_z: shifts the density function in the code
c cpumess: WALL CLOCK time required for message passing
c cpuinou: WALL CLOCK time required for I/O
c cpucomp: WALL CLOCK time required for computation
c cpus: WALL CLOCK time at start time
c cpus: WALL CLOCK time at end time


      nstart=0
      shift_c=0
      shift_l=0
      shift_z=0.0
      cpumess=0.0
      cpucomp=0.0
      cpuinou=0.0
      cpus=0.0
      cpuf=0.0


      r1n=max(i1n,r1n)
      r1x=min(i1x,r1x)
      r2n=max(i2n,r2n)
      r2x=min(i2x,r2x)
      r3n=max(i3n,r3n)
      r3x=min(i3x,r3x)

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

      end subroutine INIT_param_fortran



c rd1=3: excess data points on local node
c rd2=3: excess data points on local node
c rd3=3: excess data points on local node

c lengthx: x extension of simulation box in m
c lengthy: y extension of simulation box in m
c lengthz: z extension of simulation box in m

c i1tot: global size of the grid in x
c i2tot: global size of the grid in y
c i3tot: global size of the grid in z

c xnpe: number of PEs in x-direction
c ynpe: number of PEs in y-direction
c znpe: number of PEs in z-direction

c i1n: smallest grid index of the grid in x
c i1x: largest grid index of the grid in x
c i2n: smallest grid index of the grid in y
c i2x: largest grid index of the grid in y
c i3n: smallest grid index of the grid in z
c i3x: largest grid index of the grid in z

c i1n=i1x: 2D run in the yz-plane
c i2n=i2x: 2D run in the xz-plane
c i3n=i3x: 2D run in the xy-plane

c i1n=i1x and i2n=i2x: 1D run along z-axis
c i1n=i1x and i3n=i3x: 1D run along y-axis
c i2n=i2x and i3n=i3x: 1D run along x-axis

      subroutine INIT_param_domain

      use VLA_variables
      use PIC_variables

      implicit none

      rd1=3
      rd2=3
      rd3=3

      lengthx=1.0*1.0e-6
      lengthy=1.0*1.0e-6
      lengthz=20.0*1.0e-6

      i1tot=10.0
      i2tot=10.0
      i3tot=400.0

      xnpe=1
      ynpe=1
      znpe=1

      i1n=8
      i1x=8
      i2n=8
      i2x=8
      i3n=0
      i3x=399

      boundary_part_x=1
      boundary_part_y=1
      boundary_part_z=1

! PML BOUNDARY CONDITIONS 

! condition = 0 : radiating boundaries
! condition = 1 : periodic boundaries
! condition = 2 : upmls
! condition = 3 : radiating boundaries with upmls at a given time

      boundary_fld_x1 = 1
      boundary_fld_x2 = 1
      boundary_fld_y1 = 1
      boundary_fld_y2 = 1
      boundary_fld_z1 = 1
      boundary_fld_z2 = 1

! thick: thickness of pml in grid points
! cushion: thickness of buffer region in gridpoints
! size: thickness of pml and buffer region in gridpoints
! pml: polynomial order of pml

      thick = 10
      cushion = int(thick/3)
      size = thick+cushion
      pml = 3
       
      end subroutine INIT_param_domain

c     initialize some parameters

      subroutine INIT_param_psc

      use VLA_variables
      use PIC_variables

      implicit none

c qq: reference charge in As
c mm: reference mass in kg
c tt: reference temperature in J
c c: light velocity in m/s
c eps0: eps0 in As/Vm

c nmax: maximum number of time steps
c cpum: CPU time limit
c lw: laser wavelength in m
c i0: irradiance in W um**2/m**2
c n0: peak atom density in m**3

      qq=1.6021e-19
      mm=9.1091e-31
      tt=1.6021e-16
      cc=3.0e8
      eps0=8.8542e-12

      nmax=0 ! calculate later
      cpum=100000.0
      lw=3.2*1.0e-6
      i0=1.0e21          ! modified from 22
      n0=1.0e26          ! modified from 28
      e0=0 ! calculate later

      nicell=200         ! maximum number of particles per spezies and cell

      end subroutine INIT_param_psc

c     wrapper for C to be access these parameters

      subroutine GET_param_domain(length, itot, in, ix,
     &     bnd_fld_lo, bnd_fld_hi, bnd_part, nproc, rd)

      use VLA_variables
      use PIC_variables

      implicit none

      real(kind=8) :: length(3)
      integer :: itot(3), in(3), ix(3)
      integer :: bnd_fld_lo(3), bnd_fld_hi(3), bnd_part(3)
      integer :: nproc(3), rd(3)

      length(1)=lengthx
      length(2)=lengthy
      length(3)=lengthz

      itot(1)=i1tot
      itot(2)=i2tot
      itot(3)=i3tot

      in(1)=i1n
      in(2)=i2n
      in(3)=i3n

      ix(1)=i1x
      ix(2)=i2x
      ix(3)=i3x

      nproc(1)=xnpe
      nproc(2)=ynpe
      nproc(3)=znpe

      bnd_fld_lo(1)=boundary_fld_x1
      bnd_fld_lo(2)=boundary_fld_y1
      bnd_fld_lo(3)=boundary_fld_z1
      bnd_fld_hi(1)=boundary_fld_x2
      bnd_fld_hi(2)=boundary_fld_y2
      bnd_fld_hi(3)=boundary_fld_z2

      bnd_part(1)=boundary_part_x
      bnd_part(2)=boundary_part_y
      bnd_part(3)=boundary_part_z

      rd(1)=rd1
      rd(2)=rd2
      rd(3)=rd3

      end subroutine GET_param_domain

c     wrapper for C to be able to set these parameters

      subroutine SET_param_domain(length, itot, in, ix,
     &     bnd_fld_lo, bnd_fld_hi, bnd_part, nproc, rd)

      use VLA_variables
      use PIC_variables

      implicit none

      real(kind=8) :: length(3)
      integer :: itot(3), in(3), ix(3)
      integer :: bnd_fld_lo(3), bnd_fld_hi(3), bnd_part(3)
      integer :: nproc(3), rd(3)

      lengthx=length(1)
      lengthy=length(2)
      lengthz=length(3)

      i1tot=itot(1)
      i2tot=itot(2)
      i3tot=itot(3)

      i1n=in(1)
      i2n=in(2)
      i3n=in(3)

      i1x=ix(1)
      i2x=ix(2)
      i3x=ix(3)

      xnpe=nproc(1)
      ynpe=nproc(2)
      znpe=nproc(3)

      boundary_fld_x1=bnd_fld_lo(1)
      boundary_fld_y1=bnd_fld_lo(2)
      boundary_fld_z1=bnd_fld_lo(3)
      boundary_fld_x2=bnd_fld_hi(1)
      boundary_fld_y2=bnd_fld_hi(2)
      boundary_fld_z2=bnd_fld_hi(3)

      boundary_part_x=bnd_part(1)
      boundary_part_y=bnd_part(2)
      boundary_part_z=bnd_part(3)

      rd1=rd(1)
      rd2=rd(2)
      rd3=rd(3)

      end subroutine SET_param_domain

c     wrapper for C to be able to access these parameters

      subroutine GET_param_psc(qq_, mm_, tt_, cc_, eps0_,
     *     nmax_, cpum_, lw_, i0_, n0_, e0_, b0_, j0_, rho0_, phi0_,
     *     a0_, nicell_)

      use VLA_variables
      use PIC_variables

      implicit none

      real(kind=8) :: qq_, mm_, tt_, cc_, eps0_
      real(kind=8) :: cpum_, lw_
      real(kind=8) :: i0_, n0_, e0_, b0_, j0_, rho0_, phi0_, a0_
      integer :: nmax_, nicell_

      qq_ = qq
      mm_ = mm
      tt_ = tt
      cc_ = cc
      eps0_ = eps0

      nmax_ = nmax
      cpum_ = cpum
      lw_ = lw
      i0_ = i0
      n0_ = n0
      e0_ = e0
      b0_ = b0
      j0_ = j0
      rho0_ = rho0
      phi0_ = phi0
      a0_ = a0
      nicell_ = nicell

      end subroutine GET_param_psc

c     wrapper for C to be able to set these parameters

      subroutine SET_param_psc(qq_, mm_, tt_, cc_, eps0_,
     *     nmax_, cpum_, lw_, i0_, n0_, e0_, b0_, j0_, rho0_, phi0_,
     *     a0_, nicell_)

      use VLA_variables
      use PIC_variables

      implicit none

      real(kind=8) :: qq_, mm_, tt_, cc_, eps0_
      real(kind=8) :: cpum_, lw_
      real(kind=8) :: i0_, n0_, e0_, b0_, j0_, rho0_, phi0_, a0_
      integer :: nmax_, nicell_

      qq = qq_
      mm = mm_
      tt = tt_
      cc = cc_
      eps0 = eps0_

      nmax = nmax_
      cpum = cpum_
      lw = lw_
      i0 = i0_
      n0 = n0_
      e0 = e0_
      b0 = b0_
      j0 = j0_
      rho0 = rho0_
      phi0 = phi0_
      a0 = a0_
      nicell = nicell_

      end subroutine SET_param_psc

c     wrapper for C to be able to set these coefficients

      subroutine SET_param_coeff(cori_, alpha_, beta_, eta_,
     *     wl_, ld_, vos_, vt_, wp_, dx_, dt_, np_, nnp_)

      use VLA_variables
      use PIC_variables

      implicit none

      real(kind=8) :: cori_, alpha_, beta_, eta_
      real(kind=8) :: wl_, ld_, vos_, vt_, wp_, dt_
      real(kind=8) :: dx_(3)
      integer :: np_, nnp_

      cori = cori_
      alpha = alpha_
      beta = beta_
      eta = eta_
      wl = wl_
      ld = ld_
      vos = vos_
      vt = vt_
      wp = wp_

      dx = dx_(1)
      dy = dx_(2)
      dz = dx_(3)
      dt = dt_

      np = np_
      nnp = nnp_

      end subroutine SET_param_coeff
