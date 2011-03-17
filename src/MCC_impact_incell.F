c PERFORM IMPACT IONIZATION EVENTS IN EACH CELL

      subroutine MCC_impact_incell

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
