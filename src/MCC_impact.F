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
               if(mnq.EQ.mpart(mat)) goto 122                    ! DETERMINE MATERIAL VIA MASS
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
c            open(11,file=trim(data_out)
c     &           //'/'//node//'impact'//label,
c     &           access='sequential',form='formatted')

c            write(11,*) 'IMPACT PARAMETERS AT TIMESTEP ',n
c            write(11,*) 'unsrt p:   = ', n_unsorted
c            write(11,*) '<p_0>      = ', p_0_sum/p_0_count
c            write(11,*) 'ptoo_large = ', p_0_toolarge
            
c            close(11)

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
            if(no_ionization_yet.EQ.1.AND.mcc_np(mat).gt.0) then
               
               do ics=0, max_ics-1                             ! ics  = random charge state index

                  cs=mcc_cslist(ics)                           ! cs   = charge state
               if(no_ionization_yet.EQ.1.AND.mcc_nc(mat,cs).gt.0) then

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
                        
                        if (lph_n.gt.nialloc) then
                           write(*,*) 'ENLARGE ARRAY====='
                           call SERV_systime(cpuc)
                           call SERV_labelgen(mpe,node)
                           write(6,*) node                 
                           open(11,file=trim(data_out)//'/'
     &                          //node//'ENLARGE', 
     &                          access='sequential',form='unformatted')
                           do k=0,11*lph+10,100
                              le=min(k+99,11*lph+10)
                              write(11) (p_niloc(la),la=k,le)
                           enddo
                           close(11)
                           
                           nialloc=int(1.2*lph_n+12)
                           deallocate(p_niloc)
                           allocate(p_niloc(0:11*nialloc+10))
                           
                           open(11,file=trim(data_out)//'/'
     &                          //node//'ENLARGE',
     &                          access='sequential',form='unformatted')
                           do k=0,11*lph+10,100
                              le=min(k+99,11*lph+10)
                              read(11) (p_niloc(la),la=k,le)
                           enddo
                           close(11)
                           call SERV_systime(cpud)
                           s_cpud=s_cpud+cpud-cpuc
                           
                        endif

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





