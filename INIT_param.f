c PARAMETER INITIALIZATION 

      subroutine INIT_param

      use VLA_variables
      use PIC_variables

      implicit none


c PARAMETERS THAT SHOULD NOT BE CHANGED UNDER NORMAL CONDITIONS

c rd1=3: excess data points on local node
c rd2=3: excess data points on local node
c rd3=3: excess data points on local node
c qq: reference charge in As
c mm: reference mass in kg
c tt: reference temperature in J
c c: light velocity in m/s
c eps0: eps0 in As/Vm
c pi: pi


      rd1=3
      rd2=3
      rd3=3
      qq=1.6021e-19
      mm=9.1091e-31
      tt=1.6021e-16
      cc=3.0e8
      eps0=8.8542e-12
      pi=3.1415927


c THESE PARAMETRES ARE USED FOR DETERMINING GRID SIZE,
c PARTICLE NUMBER, DATA PARTITIONING, AND PHYSICAL
c PARAMETERS!

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

c boundary_field_x: 0=radiating, 1=periodic in x
c boundary_field_y: 0=radiating, 1=periodic in y
c boundary_field_z: 0=radiating, 1=periodic in z
c boundary_part_x: 0=reflecting, 1=periodic in x
c boundary_part_y: 0=reflecting, 1=periodic in y
c boundary_part_z: 0=reflecting, 1=periodic in z

c nmax: maximum number of time steps
c cpum: CPU time limit
c lw: laser wavelength in m
c i0: irradiance in W um**2/m**2
c n0: peak atom density in m**3


      lengthx=0.1*1.0e-6
      lengthy=0.1*1.0e-6
      lengthz=0.1*1.0e-6

      i1tot=10.0
      i2tot=10.0
      i3tot=10.0

      xnpe=1
      ynpe=1
      znpe=1

      i1n=4
      i1x=4
      i2n=0
      i2x=9
      i3n=0
      i3x=9

      boundary_field_x=1
      boundary_field_y=1
      boundary_field_z=1
      boundary_part_x=1
      boundary_part_y=1
      boundary_part_z=1

      nmax=10
      cpum=100000.0
      lw=1.0*1.0e-6
      i0=1.0e22          ! modified from 22
      n0=1.0e34          ! modified from 28


! PML BOUNDARY CONDITIONS 

! condition = 'true' : derive pml coefficients 
! condition = 'false' : do nothing
! condition = 'time' : derive pml after estimated time

      boundary_pml_x1 = 'false'
      boundary_pml_x2 = 'false'
      boundary_pml_y1 = 'false'
      boundary_pml_y2 = 'false'
      boundary_pml_z1 = 'false'
      boundary_pml_z2 = 'false'

! thick: thickness of pml in grid points
! cushion: thickness of buffer region in gridpoints
! size: thickness of pml and buffer region in gridpoints
! pml: polynomial order of pml

      thick = 10
      cushion = int(thick/3)
      size = thick+cushion
      pml = 3
       

c REQUIRED SETTINGS
c Please do not change

! checking condition for boundary setting of em fields
! no pml for periodic condition

      if (boundary_field_x==1) then
         boundary_pml_x1 = 'false'
         boundary_pml_x2 = 'false'
      end if
      if (boundary_field_y==1) then
         boundary_pml_y1 = 'false'
         boundary_pml_y2 = 'false'
      end if
      if (boundary_field_z==1) then
         boundary_pml_z1 = 'false'
         boundary_pml_z2 = 'false'
      end if


      if (i1x==i1n) then
         xnpe=1
         boundary_part_x=1
         boundary_field_x=1
         boundary_pml_x1 = 'false'
         boundary_pml_x2 = 'false'
      else
         if (i1x-i1n+1.le.(rd1+1)*xnpe) then
            i1x=i1n+(rd1+1)*xnpe-1
         endif
         if (i1x-i1n+1.le.2*(size+1)) then     ! min total size: 2*(size+1)
            i1x=i1n+2*(size+1)-1
         endif
         if (i1x-i1n+1.le.(size+1)*xnpe) then  ! min size of one node
            i1x=i1n+(size+1)*xnpe-1
         endif
         if (i1x-i1n+1.gt.i1tot) then
            i1tot=i1x-i1n+1
         endif
      endif
      if (i2x==i2n) then
         ynpe=1
         boundary_part_y=1
         boundary_field_y=1
         boundary_pml_y1 = 'false'
         boundary_pml_y2 = 'false' 
      else
         if (i2x-i2n+1.le.(rd2+1)*ynpe) then
            i2x=i2n+(rd2+1)*ynpe-1
         endif
         if (i2x-i2n+1.le.2*(size+1)) then  
            i2x=i2n+2*(size+1)-1
         endif
         if (i2x-i2n+1.le.(size+1)*ynpe) then
            i2x=i2n+(size+1)*ynpe-1
         endif
         if (i2x-i2n+1.gt.i2tot) then
            i2tot=i2x-i2n+1
         endif
      endif
      if (i3x==i3n) then
         znpe=1
         boundary_part_z=1
         boundary_field_z=1
         boundary_pml_z1 = 'false'
         boundary_pml_z2 = 'false'
      else
         if (i3x-i3n+1.le.(rd3+1)*znpe) then
            i3x=i3n+(rd3+1)*znpe-1
         endif
         if (i3x-i3n+1.le.2*(size+1)) then
            i3x=i3n+2*(size+1)-1
         endif
         if (i3x-i3n+1.le.(size+1)*znpe) then
            i3x=i3n+(size+1)*znpe-1
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
      e0=dsqrt(2.0*i0/eps0/cc)/lw/1.0e6
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

      dx=lengthx/(i1tot-1)/ld
      dy=lengthy/(i2tot-1)/ld
      dz=lengthz/(i3tot-1)/ld

      dt=0.75*sqrt(1.0/(1.0/dx**2+1.0/dy**2+1.0/dz**2))
      nnp=int(2.0*pi/dt)+1
      dt=2.0*pi/nnp
      np=2*nnp


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
      dnprf=100                                      ! changed by ab
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
      dnprparti=1
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

c      if (mpe.eq.0) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.1) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.2) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.3) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.4) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.5) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.6) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.7) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.8) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.9) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.10) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.11) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.12) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.13) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.14) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.15) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.16) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.17) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.18) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.19) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.20) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.21) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.22) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.23) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.24) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.25) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.26) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.27) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.28) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.29) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.30) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif
c      if (mpe.eq.31) then
c         data_out='/scratch/local/artb/test'
c         data_chk='/scratch/local/artb/test'
c      endif


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
      write(6,*) 'PARAM',i1x,i2x,i3x
      end subroutine INIT_param
