c ====================================================================
c ELECTRON IMPACT IONIZATION OF IONS by Andreas Kemp and Hartmut Ruhl 10/2003
c 
c Implementation by A.Kemp 03/11/2005 last modified 08/04/2005 by ak
c
c Initialize this routine by calling INIT_MCC
c ====================================================================




c     Modified to get result in SI units by Toma Toncian and Monika Omieczynski 
c     for MCC_impact .f
c     15.12.2006
c###################################################################
c
c      subroutine albeli(pe, pcf, kncf, pxs, kermsg)
       subroutine albeli(ape, pcf, kncf, pxs)
c
c     electron impact ionization cross section fits
c
c     reference: k. l. bell et al, j. phys. chem. ref. data 12, 891
c                (1983)
c
c     this is an iaea subroutine to calculate cross sections for
c     projectile energy (ev).
c
c     pe =  electron energy (ev)
c
c     the number of fitting parameters varies depending on the
c     number of terms taken in the numerical fitting and on the
c     allowance for excitation autoionization in the cross section
c     to fit cross sections with excitation autoionization two seperate
c     fits are defined. one from the ionization threshold and a second
c     fit for energies above the autoionization threshold.
c     the number of parameters in any entry is given by kncf
c
c     pcf(1) = ionization potential (ev)
c     pcf(2-7) = fitting parameters ( can be less than 6 parameters)
c
c     if cross section has excitation autoionization structure then
c     for the second fit
c
c     pcf(8) = autoionization threshold (ev)
c     pcf(9) = ionization potential (ev)
c     pcf(10-15) = fitting parameters for this fit (can be less than
c                  6 parameters)
c
c     kncf = number of parameters supplied in pcf (must be 8)
c     pxs = ionization cross section (m[2])
c     kermsg = error message, ' ' is ok
c
c     written by j. j. smith , iaea atomic and molecular data unit
c     
c
c======================================================================
c
      double precision ape, pcf, pxs
      double precision ion, power, power1, xs, a, x, x2
      dimension pcf(15)
c      character*(*) kermsg
c
c      kermsg = ' '
      if(ape .lt. pcf(1) .and. ape .eq. 0) then
        pxs = 0.0d0
        return
c     else
c        kermsg = ' '
      endif
c
c---  determine parameters to be used
c
      if (kncf .gt. 7 .and. pe .gt. pcf(8) ) then
c
c---      autoionization included and energy > autoionization threshold
c
        ion = pcf(9)
        a=pcf(10)
        istart=11
        iend=kncf
      else
        ion = pcf(1)
        a=pcf(2)
        istart=3
        if (kncf .gt. 7) then
          iend=7
        else
          iend = kncf
        endif
      endif
c
c---  generate cross section
c
      x=ion/ape
      x2= 1.0d0/x
c
c---  contribution from bethe term
c
      xs = a*dlog(x2)
      if ( kncf .ge. istart) then
c
c---  contribution from least squares fit terms
c
        power1 = 1.0d0 - x
        power = power1
        do 10 i=istart,iend
          xs = xs + pcf(i)*power
          power = power*power1
  10    continue
      endif
c
c---  scale results to m[2]
c
      pxs= 1.0d-17*xs /(ape*ion)
c
c      write(*,*) 'energie', ape, 'sigma', pxs
      return
      end subroutine albeli

c ====================================================================
c ELECTRON IMPACT IONIZATION OF IONS by Andreas Kemp and Hartmut Ruhl 10/2003
c 
c Implementation by A.Kemp 03/11/2005 last modified 08/04/2005 by ak
c
c Initialize this routine by calling INIT_MCC
c  
c To be executed after sorting particles and binary collisions,  
c    since the routine adds particles and thus messes up particle order 
c    which is required for the binary collisions
c The sort routine sometimes gets one particle in the wrong position (as of 03/2005)
c    according to the cell ordering. This leads to a small error which we accept.
c We assume that the ions are at rest in Lab frame
c
c
c
c TO DO: 1. randomize charge states and materials: in this version there is a preference to 
c           ionize low charge states first
c        2. cpuserv times off 
c        3. introduce Lotz formula for impact x-sections
c
c ====================================================================

c MAIN ROUTINE FOR ELECTRON IMPACT IONIZATION

      subroutine MCC_impact

      use PIC_variables
      use VLA_variables
      use MCC_variables

      implicit none

      character*(5) :: node, label

      integer :: qnqi
      integer :: nel,nion,nc,lnh1,lnh2
      integer :: mat,cs,m
      integer :: lnh

      real(kind=8) :: qnq,mnq
      real(kind=8) :: s_cpub
      real(kind=8) :: s_cpud
      real(kind=8) :: s_cpuf


      cpua=0.0
      cpub=0.0
      cpuc=0.0
      cpud=0.0
      cpue=0.0
      cpuf=0.0


      s_cpub=0.0
      s_cpud=0.0
      s_cpuf=0.0

      call SERV_systime(cpue)

c BEGIN LOOP OVER PARTICLES 

! note that the species index sp stands for a distinct mass and charge state, 
! not for a particular material; the meaning of sp changes in each cell. 
! sp is to distinguish the species for binary collisions,
! while the material index mat distinguishes material types for atomic physics

      if(niloc.gt.1) then 

         mcc_matlist=0          ! randomized lists of material and charge state
         mcc_cslist =-1         ! charge state zero means neutral particle
         max_imat   =0          ! maximum index value of mcc_matlist
         max_ics    =0          ! maximum index value of mcc_cslist

         mcc_elist=0
         mcc_ilist=0
         mcc_np  =0
         mcc_nc  =0             ! number of particles of material mat and charge state qnq   

         n_unsorted=0           ! number of unsorted particles in array
                                ! (error in PIC_sort) for statistics
         
         nc        =0           ! total number of charged particles
         nel       =0           ! total number of electrons
         nion      =0           ! total number of ions

         p_0_count =0           ! reset
         p_0_sum   =0.0         ! average collision probability



         do l=niloc,1,-1


            qnq=p_niloc(11*l+6)
            mnq=p_niloc(11*l+7)

            if(qnq*qnq.GT.0.0) nc=nc+1                            ! number of charged particles


c DISTINGUISH ELECTRONS AND MATERIALS FOR ATOMIC PHYSICS
            
            
            do mat=0,NMAT
               if(mpart(mat)-0.1.LE.mnq.and.mnq.LE.mpart(mat)+0.1) then
                  goto 122                                               ! DETERMINE MATERIAL VIA MASS
               endif
            enddo

            write(*,*) 'MCC_impact: unknown material'
            stop

 122        continue
             
            if(mat==0) then 
               mcc_elist(nel)=l
               nel=nel+1                                         ! electron
            endif

            qnqi=int(qnq)
            if((qnqi.GE.0).AND.(qnqi.LT.n_xstable(mat))) then    ! only ions AND exclude highest ion charge state

                nion=nion+1

                do m=0,max_imat-1                                ! check new material into mcc_matlist
                   if(mat.EQ.mcc_matlist(m)) goto 124
                enddo
                mcc_matlist(max_imat)=mat                                      
                max_imat=max_imat+1

 124            continue

                do cs=0,max_ics-1                                ! check new charge state into mcc_cslist
                   if(qnqi.EQ.mcc_cslist(cs)) goto 126
                enddo
                mcc_cslist(max_ics)=qnqi                                      
                max_ics=max_ics+1

 126            continue

            	mcc_ilist(mat,qnqi,mcc_nc(mat,qnqi))=l           ! check ion into mcc_ilist
            	mcc_np(mat)=mcc_np(mat)+1
            	mcc_nc(mat,qnqi)=mcc_nc(mat,qnqi)+1

            endif

            lnh1=p_niloc(11*l+8)                                 ! cell number of particle
            lnh2=p_niloc(11*(l-1)+8)                             ! cell number next particle (count backwards)           


            if (lnh2.ne.lnh1) then                               ! no further particles in cell

               mcc_np(0)=nel                                     ! store number of electrons
               
               if (nc.gt.0) then                                 ! AT LEAST ONE CHARGED PARTICLE FOR IMPACT IONIZATION
               
                  lnh=lnh1
                  if( lnh1<lnh2) then                            ! count errors in sort routine
                     n_unsorted=n_unsorted+1
                  endif
                  
                  if(nel.gt.0.AND.nion.gt.0) then                ! ONE ELECTRON AND ONE ION AT LEAST 
c		  write(*,*) 'MCC_impact: go'
                     call MCC_impact_incell
                  endif 
                     
               endif       
      
c     RESET CELL CHARACTERIZATION AFTER PRESENT CELL IS TREATED

               mcc_matlist=0
               mcc_cslist =-1
               max_imat   =0
               max_ics    =0  

               mcc_elist=0
               mcc_ilist=0
               mcc_np   =0
               mcc_nc   =0

               nc       =0
               nel      =0
               nion     =0

            endif
         enddo                   

c END OF LOOP OVER PARTICLES


c OUTPUT OF CRITICAL IMPACT IONIZATION PARAMETERS

            call SERV_labelgen(n,label)
            call SERV_labelgen(mpe,node)
            open(11,file=trim(data_out)
     &           //'/'//node//'impact'//label,
     &           access='sequential',form='formatted')

            write(11,*) 'IMPACT PARAMETERS AT TIMESTEP ',n
            write(11,*) 'unsrt p:   = ', n_unsorted
            write(11,*) '<p_0>      = ', p_0_sum/p_0_count
            write(11,*) 'ptoo_large = ', p_0_toolarge
            
            close(11)

            call SERV_systime(cpud)
            s_cpud=s_cpud+cpud-cpuc

      endif

      end subroutine MCC_impact





c PERFORM IMPACT IONIZATION EVENTS IN EACH CELL



      subroutine MCC_impact_incell

      use PIC_variables
      use VLA_variables
      use MCC_variables

      implicit none

      integer :: m,p,nel_now,lph,lph_n,le,la
      integer :: no_ionization_yet,rix
      integer :: imat,mat,ics,cs
      character*(5) :: node

      real(kind=8) :: nu,sv,nu_0,nu_tot,p_0
      real(kind=8) :: px,py,pz,pa,et,ek
      real(kind=8) :: R,sigma,pel

      real(kind=8) :: s_cpud

c DETERMINE NULL COLLISION FREQUENCY    

      nu_0=0.0
      do imat=0, max_imat-1                       ! only materials that are present in cell

         mat=mcc_matlist(imat)
         if(mcc_np(mat).gt.0) then

            do ics=0,max_ics-1                    ! only charge states that are present in cell

               cs=mcc_cslist(ics)                 ! charge state ic

               np=mcc_nc(mat,cs)                  ! number of particles 
               nu_0=nu_0+cori*n0*np*max_sigmav(mat,cs)  
                                                  ! max_sigmav: tabulated max impact ionization 
                                                  ! cross section times velocity for material "mat" 
                                                  ! and charge state "cs" 
            enddo
         endif
      enddo
      p_0=1.0-exp(-nu_0*dtsi)                     ! IONIZATION PROBABILITY PER TIME STEP, 
                                                  ! dtsi: time step in SI units

c GO THROUGH LIST OF ALL ELECTRONS THAT CAN TAKE PART IN IMPACT IONIZATION

      p_0_sum  =p_0_sum+p_0
      p_0_count=p_0_count+1
      if(p_0>0.095) then
         p_0_toolarge=p_0_toolarge+1
      endif
               
      nel_pro=nel_pro+p_0*mcc_np(0)               ! p_0*mcc_np(0): number of ionizing electrons per cell, 
                                                  ! accumulate small probabilities
      nel_now=int(nel_pro)
      nel_pro=nel_pro-nel_now
      
      do j=0,nel_now-1                            ! for the first nel_now electrons
                                                  ! note: electrons are assumed to 
                                                  ! be randomized here
         
         p=mcc_elist(j)                           ! get electron particle index
         no_ionization_yet=1                      ! permit only one ionization event per electron
         
         call random_number(R)
         R=R*nu_0
         nu_tot=0.0
         
         do imat=0, max_imat-1                                 ! imat = random material index

            mat=mcc_matlist(imat)                              ! mat  = material number
            if(no_ionization_yet.AND.mcc_np(mat).gt.0) then
               
               do ics=0, max_ics-1                             ! ics  = random charge state index

                  cs=mcc_cslist(ics)                           ! cs   = charge state
                  if(no_ionization_yet.AND.mcc_nc(mat,cs).gt.0) then

                     px=p_niloc(11*p+3)                        ! assign e-impact momentum
                     py=p_niloc(11*p+4)                        ! assuming that ion is at rest
                     pz=p_niloc(11*p+5)
                     pa=sqrt(px*px+py*py+pz*pz)
                     et=me*sqrt(1.0+pa*pa)
                     ek=et-me

                     call MCC_ixsection(mat,cs,ek,sigma,sv)

                                                                ! ek : kinetic energy of ionizing electron, 
                                                                ! sv : sigmav, get cross section
                     nu=mcc_nc(mat,cs)*n0*cori*sv
                     nu_tot=nu_tot+nu
                     
                     if(R <= nu_tot) then                       ! PROCEED IMPACT IONIZATION EVENT:
                                              
                        m=0
                        rix=-1
                        do while (m.EQ.0)
                          rix=rix+1                             ! ions are already randomized
                          m=mcc_ilist(mat,cs,rix)               ! m=random ion particle index
                        enddo

                        mcc_ilist(mat,cs,rix)=0                 ! exclude ionzed particle from list
                        mcc_nc(mat,cs)=mcc_nc(mat,cs)-1         ! update lists

                        p_niloc(11*m+6)=p_niloc(11*m+6)+1.0     ! increase ion charge state
                        
                        px=px/pa
                        py=py/pa                                ! normalize electron momentum
                        pz=pz/pa 
                        
                        et=et-xstable_t(mat,cs)              ! subtract ionization energy
                        pel=sqrt(et*et-me*me)/me               ! from electron kinetic energy..
                        
                        px=px*pel                              ! ..and reduce its momentum 
                        py=py*pel                              ! to account for 
                        pz=pz*pel                              ! the ionization enery
                                                               ! NOTE: no momentum conservation but 
                                                               ! negligible for small i-poten's

                        lph=niloc                              ! create a new electron:
                        lph_n=lph+1

                        call REALLOC_particles(lph_n)

                        lph=lph_n
                        niloc=lph
                        
                        p_niloc(11*lph+0)=p_niloc(11*m+0)        ! assign position of 
                        p_niloc(11*lph+1)=p_niloc(11*m+1)        ! mother ion to new-born electron
                        p_niloc(11*lph+2)=p_niloc(11*m+2)
                        p_niloc(11*lph+3)=p_niloc(11*m+3)        ! initialize with ion velocity
                        p_niloc(11*lph+4)=p_niloc(11*m+4)
                        p_niloc(11*lph+5)=p_niloc(11*m+5)
                        p_niloc(11*lph+6)=-1.0d0
                        p_niloc(11*lph+7)=+1.0d0
                        p_niloc(11*lph+8)=p_niloc(11*m+8)        ! assign local cell number
                        p_niloc(11*lph+9)=p_niloc(11*m+9)        ! 
                        p_niloc(11*lph+10)=p_niloc(11*m+10)      ! 

                        no_ionization_yet=0                      ! allow only one ionization event per el
                     endif

                  endif
               enddo

            endif
         enddo
      enddo

      end subroutine MCC_impact_incell
      


c CALCULATE IONIZATION X-SECTION FROM LOTZ FORMULA TO BE SIMPLE AND GENERAL

      subroutine MCC_ixsection(mat,xs,ek,sigma,sv)! interpolate x-section and sv == sigma v
                                                  ! for material mat and initial charge state xs
                                                  ! expect ek in eV
                                                  ! return sigma and sv in SI units
      use PIC_variables
      use VLA_variables
      use MCC_variables

      implicit none

      integer      :: mat,xs                      ! expect real material and charge state index -- NOT mcc_list - index
      real(kind=8) :: sigma,sv,ek

      real(kind=8) :: va
      integer      :: ux,lx,mx

      va=cc*sqrt(ek*ek+2.0*me*ek)/(ek+me)         ! va given in SI units
      ux=xstable_n(mat,xs)-1                      ! assume table is ordered
      lx=0                                        ! for increasing energies

      if (ek < xstable_e(mat,xs,lx)) then         ! ekin < table boundary
         sigma=0
         sv=0
         return
      endif 

      if (ek > xstable_e(mat,xs,ux)) then         ! ekin > table boundary
         sigma=xstable_s(mat,xs,ux)
         sv=sigma*va
         return
      else 
         do while (ux > lx+1)                     ! perform linear bisection to find sigma
            mx=int((lx+ux)/2.0)
            if (xstable_e(mat,xs,mx) > ek) then
               ux=mx
            else
               lx=mx
            endif
         enddo
         sigma=((ek-xstable_e(mat,xs,lx))*xstable_s(mat,xs,ux)+
     c        (xstable_e(mat,xs,ux)-ek)*xstable_s(mat,xs,lx))/
     c        (xstable_e(mat,xs,ux)-xstable_e(mat,xs,lx))
         sv=sigma*va
         return
      endif

      end subroutine MCC_ixsection




c     WRITE IONIZATION X-SECTIONS INTO A FILE FOR CONTROL PURPOSES ONLY

      subroutine MCC_write_ixsections

      use MCC_variables
      
      implicit none

      real(kind=8) :: e0,e1,f,ek
      real(kind=8) :: sigma,sv
      integer      :: i,j,mat
      integer      :: NPOINTS

      NPOINTS=100

      do mat=1,NMAT

         write(*,*) '# Electron Impact Ionization X-Section for : ',
     &        matname(mat)
         do i=0,n_xstable(mat)-1
            
            write(*,*) '# charge state ',i
            write(*,*) '# Ekin[eV]  sigma[m^2]  sigma v[m^3/s]'
            e0=xstable_t(mat,i)
            e1=xstable_e(mat,i,xstable_n(mat,i)-1)
            f=(e1/e0)**(1.0/NPOINTS)
            
            ek=e0
            do j=0,NPOINTS
               call MCC_ixsection(mat,i,ek,sigma,sv)
               write(*,*) ek,sigma, sv
               ek=ek*f
            enddo
            write(*,*) ''
         enddo
      enddo

      end subroutine MCC_write_ixsections





c     INITIALIZE IMPACT IONIZATION PACKAGE BEFORE CALLING MCC_impact:

      subroutine INIT_MCC

      use PIC_variables
      use VLA_variables
      use MCC_variables

      implicit none

      integer :: mat
      integer :: cs
      integer :: NDAT
      
      integer :: kncfmax 				!albeli
      integer :: kncf					!albeli
      double precision:: pcf, pxs				!albeli
      double precision ape
      integer :: counteralbeli
      dimension pcf(15)					!albeli
c      character*(*) kermsg
c     for albeli subroutine kncfmax=15 min=8
c     and pcf fit coeficents
      real(kind=8) :: ek,ekmax,sv,va

	
      NMAT          =4                                  ! number of materials
                                                        ! in case of NMAT>1, change table allocation
      NCS           =15                                  ! number of allowed charge states>0

      NPMAX         =3*(NCS+1)*nicell                   ! conservative estimate for max number of particles in cell
      
      kncfmax=15
   	
      
      allocate(mpart(0:NMAT))
      

                                                        ! MATERIALS ARE IDENTIFIED BY MASS

      write(*,*) 'MONTE-CARLO COLLISION PACKAGE:'
      write(*,*) 'NMAT = ',NMAT
      write(*,*) 'NCS  = ',NCS
      write(*,*) 'NMPAX= ',NPMAX
      
      dtsi          =dt/wl                              ! time step in SI units


      allocate(mcc_elist(0:NPMAX))
      allocate(mcc_ilist(1:NMAT,0:NCS-1,1:NPMAX))
      allocate(mcc_nc(1:NMAT,0:NCS-1))
      allocate(mcc_np(1:NMAT))

      allocate(mcc_matlist(0:NMAT))
      allocate(mcc_cslist(0:NCS))

      p_0_toolarge  =0                                  ! number of collision events where
                                                        ! the probability is too large
      nel_pro       =0.50d0                             ! number of ii events to proceed
c      NDAT          =125                                ! max number of impact x-section data points
      NDAT          =1000	     
      me            =5.11d5                             ! electron rest mass in eV

      allocate(matname(1:NMAT))



      allocate(n_xstable(1:NMAT)) 
      allocate(xstable_n(1:NMAT, 0:NCS-1))
      allocate(xstable_t(1:NMAT, 0:NCS-1))
      allocate(max_sigmav(1:NMAT, 0:NCS-1))
      allocate(xstable_e(1:NMAT, 0:NCS-1, 0:NDAT))
      allocate(xstable_s(1:NMAT, 0:NCS-1, 0:NDAT))
       

c     INITIALIZE IMPACT TABLES

      n_xstable =(/ 2, 7, 10, 15 /)                                    ! number of tables == charge states>0
      matname(:)=(/'HELIUM','NITROGEN', 'NEON', 'ARGON'/)                         ! names of materials
c     OF course there are more charge stages for NE and AR but we are lasy 		
c     DEFAULT VALUES FOR ELECTRONS

      mat=0
      mpart(mat)=1.0                                           ! electron mass

c     IMPACT IONIZATION OF HELIUM

      mat=1
      if(mat.gt.NMAT) then
         write(*,*) 'material not allowed'
         stop
      endif
	
c      n_xstable(mat)=2  
      mpart(mat)         =7344.0                               ! use ion mass to recognize material
      xstable_n(mat,:)   =(/ 1000, 1000 /)                        ! number of array elements < NDAT
      xstable_t(mat,:)   =(/ 24.59d0 , 55.0d0  /)              ! E_threshold in eV  
                                                               ! sigma(Ekin): [m^2] ([eV])
      
      cs=0                                               ! He, Z=0-->1
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8									! Use Fit-Parameters from http://www-amdis.iaea.org/GENIE/ fuer Wirkungsquerschnitte
      pcf(1:kncf)=(/ 2.46000d+01 , 5.72000d-01 , -3.44000d-01 , 		
     c -5.23000d-01, 3.445d+00 , -6.82100d+00 , 5.57800d+00 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=10.0d0+j*(10000.0-10.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=1                                               ! He, Z=1-->2
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 5.44000d+01 , 1.85000d-01 , 8.90000d-02 , 
     c 1.31000d-01 , 3.88000d-01 , -1.09100d+00 ,  1.35400d+00 ,
     c 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=40.0d0+j*(10000.0-40.0)/1000.0		! Table with 1000 entries for cross section
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)			! Generate automatically fit for cross sections via Genie-PRG Albeli 
      xstable_s(mat,cs,j)=pxs							! kncf=Number of fit parameters required by albeli (min = 8, max = 15)
      enddo 


c     IMPACT IONIZATION OF NITROGEN

      mat=2
      if(mat.gt.NMAT) then
         write(*,*) 'material not allowed'
         stop
      endif

      mpart(mat)  =25704.0 

c      n_xstable(mat)   =7                                    ! number of tables == charge states>0
      xstable_n(mat,:) =(/ 1000 , 1000 , 1000 , 1000 , 1000 ,
     c  1000 , 1000 /)                        ! number of array elements < NDAT
      xstable_t(mat,:) =(/ 14.5d0 , 29.6d0 , 47.5d0 , 77.5d0 ,
     c    97.9d0, 5.52d2 , 6.67d2 /)                           ! E_threshold in eV  
                                                              ! sigma(Ekin): [m^2] ([eV])

      cs=0                                                   !  N, Z=0-->1
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      kncf=8
      pcf(1:kncf)=(/ 1.45d1 , 2.265d0 , -1.71d0 , -2.322d0 , 1.732d0 , 
     c 0.0d0 , 0.0d0 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=10.0d0+j*(10000.0-10.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 

      cs=1                                                              ! N, Z=1-->2
       if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      kncf=8
      pcf(1:kncf)=(/ 2.96d1 , 1.076d0 , -8.29d-1 , 8.72d-1 , -1.62d-1 , 
     c 1.533d0 , 0.0d0 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=25.0d0+j*(10000.0-25.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=2                                                              ! N, Z=2-->3
       if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      kncf=8
      pcf(1:kncf)=(/ 4.75d1 , 5.0d-1 , 2.23d-1 , 2.207d0 , -4.155d-0 , 
     c 3.769d0 , 0.0d0 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=40.0d0+j*(10000.0-40.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
       cs=3                                                              ! N, Z=3-->4
       if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      kncf=8
      pcf(1:kncf)=(/ 7.75d1 , 8.130d-1 , -7.0d-3 , -4.6d-2 , 0.0d0 , 
     c 0.0d0 , 0.0d0 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=70.0d0+j*(10000.0-70.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
       cs=4                                                              ! N, Z=4-->5
       if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      kncf=15
      pcf(1:kncf)=(/  9.79000d+01 , 2.18000d-01 , 2.38000d-01 , -2.20000d-01 ,
     c -4.46000d-01 , 2.52300d+00 , -1.90200d+00 , 4.11180d+02 ,
     c 9.79000d+01 , 8.37000d-01 ,   -2.14000d-01 , -2.53800d+00 ,
     c 7.48800d+00 , -1.10060d+01 , 5.52300d+00 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=90.0d0+j*(10000.0-90.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      
          cs=5                                                             ! N, Z=5-->6
       if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      kncf=8
      pcf(1:kncf)=(/ 5.52100d+02 , 7.96000d-01 , -5.00000d-01 , 8.84000d-01 ,
     c 0, 0, 0, 0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=540.0d0+j*(10000.0-540.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      
       cs=6                                                             ! N, Z=6-->7
       if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      kncf=8
      pcf(1:kncf)=(/ 6.67000d+02  , 4.00000d-01 , 0 , 0 ,
     c 0, 0, 0, 0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=540.0d0+j*(10000.0-540.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo

c     IMPACT IONIZATION OF NEON

      mat=3
      if(mat.gt.NMAT) then
         write(*,*) 'material not allowed'
         stop
      endif
	
c      n_xstable(mat)=10  
      mpart(mat)         =36720.0                               ! use ion mass to recognize material
      xstable_n(mat,:)   =(/ 1000, 1000 , 1000 , 1000 ,
     c 1000, 1000 ,1000, 1000 ,1000, 1000 /)                        ! number of array elements < NDAT
      xstable_t(mat,:)   =(/ 2.16d1 , 4.110d1,
     c  6.35d1, 9.25d1 , 1.262d2, 1.579d2, 2.075d2, 2.391d2 ,
     c 1.196d3, 1.3606d3  /)              ! E_threshold in eV  
                                                               ! sigma(Ekin): [m^2] ([eV])
      
      cs=0                                               ! Ne, Z=0-->1
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 2.16000d+01 , 2.19200d+00 , -4.47000d-01 ,
     c  -7.00600d+00 ,5.92700d+00, 0.0d0 , 0.0d0 ,0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=10.0d0+j*(10000.0-10.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=1                                               ! Ne, Z=1-->2
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 4.11000d+01 , 2.70500d+00 , -2.94600d+00 ,
     c 4.86200d+00 , -1.50700d+01 , 1.77800d+01 , -5.71600d+00,
     c 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=40.0d0+j*(10000.0-40.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=2                                               ! Ne, Z=2-->3 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 6.35000d+01 , 3.70100d+00 , -1.12800d+00 ,
     c  -6.34400d+00 , 4.84200d+00 , 0.0d0 , 0.0d+00, 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=60.0d0+j*(10000.0-60.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
       
      cs=3                                               ! Ne, Z=3-->4 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 9.25000d+01 , 7.85000d-01 , 1.70900d+00 ,
     c  -1.08500d+01 , 4.15000d+01 , -5.80500d+01 , 3.07200d+01,
     c 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=90.0d0+j*(10000.0-90.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 

      cs=4                                               ! Ne, Z=4-->5 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 1.26200d+02 , 1.06600d+00 , 4.42000d-01
     c , 4.75000d-01 , -2.96100d+00 , 4.47000d+00, 0.0d0,
     c 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=120.0d0+j*(10000.0-120.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=5                                               ! Ne, Z=5-->6 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 1.57900d+02 , 1.04500d+00 ,  -6.52000d-01
     c , 1.29900d+00, 0.0d+00 , 0.00d+00, 0.0d0, 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=150.0d0+j*(10000.0-150.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      
      cs=6                                               ! Ne, Z=6-->7 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  2.07500d+02 , 7.37000d-01 , 3.80000d-02
     c , 2.73000d-01 ,-3.84000d-01 , 1.33000d-01 ,
     c 0.0d0, 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=200.0d0+j*(10000.0-200.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      
      cs=7                                               ! Ne, Z=7-->8
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=14
      pcf(1:kncf)=(/  2.39100d+02 , 4.20000d-02 , -1.12000d-01
     c , 2.32500d+00 ,-1.74300d+00 , 1.50000d-02 , 0.00000d+00
     c , 8.96620d+02 ,2.39100d+02 , 6.47000d-01 , -7.45000d-01
     c , 2.76800d+00 ,-4.05900d+00 ,  1.52000d+00 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=230.0d0+j*(10000.0-230.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      cs=8                                               ! Ne, Z=8-->9
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  1.19600d+03 , 8.84000d-01 , -3.48000d-01
     c , 8.72000d-01 ,-8.44300d+00  , 1.71600d+01 ,
     c -9.17900d+00,  0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=1.10d3+j*(10000.0-1.1d3)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      cs=9                                               ! Ne, Z=9-->10
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  1.36060d+03 , 4.34000d-01 , -6.90000d-02
     c  , 9.30000d-02 , -3.90000d-02  ,0.0d0 ,0.0d0,  0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=1.30d3+j*(10000.0-1.3d3)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo      

c     IMPACT IONIZATION OF ARGON

      mat=4
      if(mat.gt.NMAT) then
         write(*,*) 'material not allowed'
         stop
      endif
	
c      n_xstable(mat)=10  
      mpart(mat)         =73440.0                               ! use ion mass to recognize material
      xstable_n(mat,:)   =(/ 1000, 1000 , 1000 , 1000 ,
     c 1000, 1000 ,1000, 1000 ,1000, 1000 , 1000 ,
     c 1000 , 1000 , 1000 , 1000 /)                        ! number of array elements < NDAT
      xstable_t(mat,:)   =(/ 1.58d1 , 2.74d1,
     c  4.07d1, 5.23d1 , 7.5d1, 9.10d1, 1.243d2 ,
     c  1.435d2, 4.225d2 , 4.787d2 ,
     c  5.39d2, 6.182d2 , 6.861d2 , 7.557d2 , 8.548d2 /)              ! E_threshold in eV  
                                                               ! sigma(Ekin): [m^2] ([eV])
      
      cs=0                                               ! Ar, Z=0-->1
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=15
      pcf(1:kncf)=(/ 1.58000d+01 , 2.53200d+00 , -2.67200d+00 ,
     c  2.54300d+00 , -7.69000d-01 , 8.00000d-03 , 
     c  6.00000d-03 , 6.19500d+01 , 1.58000d+01
     c  , 4.33700d+00 , 3.09200d+00 , -2.12530d+01 ,
     c  1.46260d+01 , 1.80000d-02 , 3.10000d-02 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=10.0d0+j*(10000.0-10.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=1                                               ! Ar, Z=1-->2
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 2.74000d+01 , 2.89600d+00 , 7.77000d-01 ,
     c  -4.44700d+00 , 2.86700d+00 , 7.20000d-02,
     c   0.0d0 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=20.0d0+j*(10000.0-20.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=2                                               ! Ar, Z=2-->3 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 4.07000d+01 , 2.08600d+00 , 1.07700d+00 ,
     c  -2.17200d+00 , 8.09000d-01 , 0.0d0 , 0.0d+00, 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=38.0d0+j*(10000.0-38.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
       
      cs=3                                               ! Ar, Z=3-->4 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 5.23000d+01 , 1.18600d+00 , -1.18000d+00 
     c   , 1.10500d+01 , -3.07900d+01 , 3.66200d+01 ,
     c   -1.54200d+01 ,   0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=50.0d0+j*(10000.0-50.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 

      cs=4                                               ! Ar, Z=4-->5 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=13
      pcf(1:kncf)=(/ 7.50000d+01 , 1.57400d+00 , 7.22000d-01 ,
     c  -2.68700d+00 , 1.85600d+00 , 0.00000d+00 , 
     c  0.00000d+00 , 2.19750d+02 ,  1.26000d+02
     c , 2.79800d+00 , 4.11400d+00 , -3.10300d+00 , 
     c   4.38000d-01 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=70.0d0+j*(10000.0-70.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=5                                               ! Ar, Z=5-->6 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=13
      pcf(1:kncf)=(/ 9.10000d+01 , 1.17000d+00 , 8.43000d-01
     c  , -2.87700d+00 , 1.95800d+00 , 0.00000d+00 ,
     c  0.00000d+00 , 2.30230d+02 , 2.00000d+02
     c  , 3.77100d+00 , 1.61630d+01  ,-3.49520d+01 ,
     c  2.08530d+01 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=85.0d0+j*(10000.0-85.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      
      cs=6                                               ! Ar, Z=6-->7 
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=15
      pcf(1:kncf)=(/  1.24300d+02 , 9.68000d-01 , 
     c  -3.06000d-01 , 2.23000d-01 , -5.00000d-03  ,
     c 1.80000d-02 , 1.70000d-02 , 2.39320d+02 ,
     c  2.20000d+02  , 3.73900d+00 , 6.81700d+00 ,
     c  -1.36650d+01 , 7.35700d+00 , -1.17000d-01
     c  , -1.39000d-01 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=120.0d0+j*(10000.0-120.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      
      cs=7                                               ! Ar, Z=7-->8
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  1.43500d+02 ,  9.07000d-01 
     c  , 1.05000d-01 , 1.47000d-01 , -7.30000d-02 
     c  , -1.00000d-03 , 0.0d+00 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=140.0d0+j*(10000.0-140.0)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      cs=8                                               ! Ar, Z=8-->9
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  4.22500d+02 , 3.51200d+00 ,
     c -3.57000d-01 , 1.15800d+00  ,
     c  -7.76000d-01 , 0.0d0 , 0.0d0 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=4.10d2+j*(10000.0-4.10d2)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo
      
      cs=9                                               ! Ar, Z=9-->10
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  4.78700d+02 , 2.79400d+00 ,
     c  4.69000d-01 , -1.29400d+01  , 2.62600d+01 
     c  , -1.34300d+01 ,0.0d0,  0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=4.70d2+j*(10000.0-4.70d2)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo      
       
      cs=10                                               ! Ar, Z=10-->11
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  5.39000d+02 ,  2.01900d+00 ,
     c  -1.32000d+00 , 1.70000d+00 ,0.0d0,  0.0d0 ,
     c  0.0d0 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=5.20d2+j*(10000.0-5.20d2)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      cs=11                                               ! Ar, Z=11-->12
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/ 6.18200d+02 , 7.85000d-01 , 
     c  1.70900d+00 , -1.08500d+01 ,
     c  4.15000d+01 , -5.80500d+01 , 3.07200d+01 
     c , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=6.00d2+j*(10000.0-6.00d2)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      
      cs=12                                               ! Ar, Z=12-->13
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  6.86100d+02  , 1.06600d+00 ,
     c  4.42000d-01 , 4.75000d-01 , -2.96100d+00 
     c , 4.47000d+00 ,0.0d0,  0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=6.80d2+j*(10000.0-6.80d2)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      
      cs=13                                               ! Ar, Z=13-->14
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  7.55700d+02 , 1.20000d+00 ,
     c  -6.52000d-01 , 1.29900d+00 ,0.0d0,  0.0d0,
     c 0.0d0 , 0.0d0 /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=7.50d2+j*(10000.0-7.50d2)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      
      cs=14                                               ! Ar, Z=14-->15
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      
      kncf=8
      pcf(1:kncf)=(/  8.54800d+02 , 6.68000d-01 ,
     c   1.97000d-01 , 1.20000d-01 ,
     c  -1.78000d-01 ,0.0d0,  0.0d0 , 0.0d0  /) 
      do j=0,xstable_n(mat,cs) -1
      pxs=0
      xstable_e(mat,cs,j)=8.45d2+j*(10000.0-8.45d2)/1000.0
      
      call albeli(xstable_e(mat,cs,j),pcf,kncf,pxs)
      xstable_s(mat,cs,j)=pxs
      enddo 
      
      
      write(*,*) '# ============================'
      write(*,*) '# INITIALIZE E-IMPACT IONIZATION'

      do mat=1,NMAT
         write(*,*) '# X-SECTIONS FOR ', matname(mat)
         do cs=0,n_xstable(mat)-1

            max_sigmav(mat,cs)=0.0
            ekmax             =0.0
            do j=0,xstable_n(mat,cs)-1
               ek=xstable_e(mat,cs,j)
               va=cc*sqrt(ek*ek+2.0*me*ek)/(ek+me)
               sv=xstable_s(mat,cs,j)*va
               if(sv > max_sigmav(mat,cs)) then
                  max_sigmav(mat,cs)=sv
                  ekmax        =ek
               endif
            enddo
            write(*,*) '# Z_0 max(sigma v)[SI] @ Ek[eV]'
            write(*,*) cs, max_sigmav(mat,cs), ekmax
         enddo
      enddo
      write(*,*) '# ============================='

      call MCC_write_ixsections

      end subroutine INIT_MCC
