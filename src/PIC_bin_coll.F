c ====================================================================
c BINARY MONTE CARLO COLLISION OPERATOR by A.Kemp and H. Ruhl, 06/2003
c THIS SUBROUTINE GROUPS PARTICLES IN PAIRS AND SCATTERS THEM.THE
c ROUTINE MUST BE SET UP PROPERLY FOR EACH SIMULATION TO TAKE ACCOUNT
c OF ALL POSSIBLE BINARY COLLISIONS.
c ====================================================================

      module PIC_bin_collision

      use PSC_domain, only: mpe

c nprc: first time step for output of collision parameters
c dnprc: timestep increment for output of collision parameters

      integer :: nprc = 0
      integer :: dnprc = 100

      real(kind=8) :: nudt,coll_add,coll_tot
      real(kind=8), dimension(:,:), allocatable :: nudt_min
      real(kind=8), dimension(:,:), allocatable :: nudt_max
      real(kind=8), dimension(:,:), allocatable :: nudt_mean
      real(kind=8), dimension(:,:), allocatable :: coll_count

      end module PIC_bin_collision
 

      subroutine PIC_bin_coll(niloc, p_niloc)

      use PIC_bin_collision
      use VLA_variables, only: n
      use OUT_params
      use PSC_globals
      use PSC_params

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)

      integer :: i,j,k,l,sp,sp1,sp2,spx_nsp,ierr
      integer :: intg,limit
      integer :: cnh1,cnh2
      integer :: spx_all,NPMAX

      character*(5) :: node,label

      real(kind=8) :: qni,mni,wni
      real(kind=8) :: rest,rho,lnL

      integer,allocatable,dimension(:,:) :: spx_list
      integer,allocatable,dimension(:) :: spx_np


      lnL=5.0d0                                   ! Coulomb logarithm
      NPMAX=5000                                  ! max number of particles per cell
      spx_nsp=7                                   ! determines maximum number of species in cell

      allocate(spx_list(0:spx_nsp,0:NPMAX))       ! matrix spx_list required for selection rules
      allocate(spx_np(0:spx_nsp))

      spx_list=0                                  ! index list
      spx_np=0                                    ! number of particles of specie sp in cell

      allocate(coll_count(0:spx_nsp,0:spx_nsp))
      allocate(nudt_min(0:spx_nsp,0:spx_nsp))
      allocate(nudt_max(0:spx_nsp,0:spx_nsp))
      allocate(nudt_mean(0:spx_nsp,0:spx_nsp))

      coll_tot=0.0d0                              ! coll total
      coll_count=0.0d0                            ! coll counter
      nudt_min=1.0d30                             ! coll freq min
      nudt_max=0.0d0                              ! coll freq max
      nudt_mean=0.0d0                             ! coll freq mean


      if (niloc.gt.1) then


c BEGIN LOOP OVER PARTICLES


         do l=niloc,1,-1

            qni=p_niloc(11*l+6)                                ! mass of particle
            mni=p_niloc(11*l+7)                                ! charge of particle

            if (qni==-1.0.and.mni==1.0) then
               sp=0                                            ! electron
               spx_list(sp,spx_np(sp))=l                       ! store particle index in list
               spx_np(sp)=spx_np(sp)+1                         ! increase particle counter by one

               if (spx_np(sp).gt.NPMAX) then
                  if (mpe==0) then
                     write(6,*) "PIC_bin_coll: NPMAX too small!"
                  endif
                  call MPI_FINALIZE(ierr)
                  stop 
               endif

            else if (qni==+1.0.and.mni==1836.0) then
               sp=1                                            ! H1+
               spx_list(sp,spx_np(sp))=l                       ! store particle index in list
               spx_np(sp)=spx_np(sp)+1                         ! increase particle counter by one

               if (spx_np(sp).gt.NPMAX) then
                  if (mpe==0) then
                     write(6,*) "PIC_bin_coll: NPMAX too small!"
                  endif
                  call MPI_FINALIZE(ierr)
                  stop 
               endif

            else if (qni==+2.0.and.mni==1836.0) then
               sp=2                                            ! H2+
               spx_list(sp,spx_np(sp))=l                       ! store particle index in list
               spx_np(sp)=spx_np(sp)+1                         ! increase particle counter by one

               if (spx_np(sp).gt.NPMAX) then
                  if (mpe==0) then
                     write(6,*) "PIC_bin_coll: NPMAX too small!"
                  endif
                  call MPI_FINALIZE(ierr)
                  stop 
               endif

            else if (qni==+1.0.and.mni==4*1836.0) then
               sp=3                                            ! He1+
               spx_list(sp,spx_np(sp))=l                       ! store particle index in list
               spx_np(sp)=spx_np(sp)+1                         ! increase particle counter by one

               if (spx_np(sp).gt.NPMAX) then
                  if (mpe==0) then
                     write(6,*) "PIC_bin_coll: NPMAX too small!"
                  endif
                  call MPI_FINALIZE(ierr)
                  stop 
               endif

            else if (qni==+2.0.and.mni==4*1836.0) then
               sp=4                                            ! He2+
               spx_list(sp,spx_np(sp))=l                       ! store particle index in list
               spx_np(sp)=spx_np(sp)+1                         ! increase particle counter by one

               if (spx_np(sp).gt.NPMAX) then
                  if (mpe==0) then
                     write(6,*) "PIC_bin_coll: NPMAX too small!"
                  endif
                  call MPI_FINALIZE(ierr)
                  stop 
               endif

            else if (qni==+1.0.and.mni==12*1836.0) then
               sp=5                                            ! C+
               spx_list(sp,spx_np(sp))=l                       ! store particle index in list
               spx_np(sp)=spx_np(sp)+1                         ! increase particle counter by one

               if (spx_np(sp).gt.NPMAX) then
                  if (mpe==0) then
                     write(6,*) "PIC_bin_coll: NPMAX too small!"
                  endif
                  call MPI_FINALIZE(ierr)
                  stop 
               endif

            else if (qni==+2.0.and.mni==12*1836.0) then
               sp=6                                            ! C2+
               spx_list(sp,spx_np(sp))=l                       ! store particle index in list
               spx_np(sp)=spx_np(sp)+1                         ! increase particle counter by one

               if (spx_np(sp).gt.NPMAX) then
                  if (mpe==0) then
                     write(6,*) "PIC_bin_coll: NPMAX too small!"
                  endif
                  call MPI_FINALIZE(ierr)
                  stop 
               endif

            else if (qni==+3.0.and.mni==12*1836.0) then
               sp=7                                            ! C3+
               spx_list(sp,spx_np(sp))=l                       ! store particle index in list
               spx_np(sp)=spx_np(sp)+1                         ! increase particle counter by one

               if (spx_np(sp).gt.NPMAX) then
                  if (mpe==0) then
                     write(6,*) "PIC_bin_coll: NPMAX too small!"
                  endif
                  call MPI_FINALIZE(ierr)
                  stop 
               endif

            endif


c***********************************************************************************************
c EXAMPLE FOR THE MATRIX spx_list 

c              sp/i   | 0  | 1  | 2  | ...| spx_nsp      sp is along x-axis, i is along y-axis, where
c           ---------------------------------------      spx_np is the number of particles of specie sp. 
c               0     | 35 | 41 | 45 | ...| 0
c               1     | 36 | 42 | 0  | ...| 0            The values in the matrix spx_list are the
c               2     | 37 | 43 | 0  | ...| 0            particle labels.
c               3     | 38 | 44 | 0  | ...| 0
c               4     | 39 | 0  | 0  | ...| 0            spx_np(0)=6
c               5     | 40 | 0  | 0  | ...| 0            spx_np(1)=4
c               6     | 0  | 0  | 0  | ...| 0            spx_np(2)=1
c              ...    | ...| ...| ...| ...| ...
c              NPMAX. | 0  | 0  | 0  | ...| 0

c***********************************************************************************************

c PAIR SELECTION RULES

c INTRA SPECIES              INTER SPECIES
c spx_np(sp1) is even        spx_np(sp1)=spx_np(sp2)
         
c            (sp1,0)         (sp1,0) - (sp2,0)
c               |            (sp1,1) - (sp2,1)
c            (sp1,1)         (sp1,2) - (sp2,2)
c--------------------
c            (sp1,2)
c               |
c            (sp1,3)
c--------------------       spx_np(sp1)>spx_np(sp2)
c            (sp1,4)
c               |           (sp1,0) - (sp2,0)    (sp1,0) - (sp2,0)    (sp1,0) - (sp2,0)    (sp1,0) - (sp2,0)
c            (sp1,5)        (sp1,1) - (sp2,0)    (sp1,1) - (sp2,0)    (sp1,1) - (sp2,0)    (sp1,1) - (sp2,0)
c                           (sp1,2) - (sp2,1)    (sp1,2) - (sp2,1)    (sp1,2) - (sp2,0)    (sp1,2) - (sp2,1)
c spx_np(sp1) is odd                             (sp1,3) - (sp2,1)    (sp1,3) - (sp2,1)    (sp1,3) - (sp2,1)
c                                                                     (sp1,4) - (sp2,1)    (sp1,4) - (sp2,2)
c            (sp1,0)
c               |
c            (sp1,1)         spx_np(sp1)<spx_np(sp2)
c--------------------
c            (sp1,0)         (sp1,0) - (sp2,0)
c               |            (sp1,0) - (sp2,1)    Here sp1 is exchanged with sp2! 
c            (sp1,2)         (sp1,1) - (sp2,2)
c--------------------
c            (sp1,1)
c               |
c            (sp1,2)         The idea is to avoid as many multiple assignments as possible. This is
c--------------------        achieved by dividing the larger particle number by the smaller one and
c            (sp1,3)         to look for the next nearest integers.that smaller and larger than this
c               |            value. Representative examples are given above. Multiple assignments do
c            (sp1,4)         not pose a problem since the Boltzmann collision operator scatters each 
c--------------------        particle with all the other particles. However, they increase the work
c            (sp1,5)         load.
c               |
c            (sp1,6)

c***********************************************************************************************



            cnh1=p_niloc(11*l+8)                                ! determine cell number of particle
            cnh2=p_niloc(11*(l-1)+8)                            ! determine cell number of next particle
            
            if (cnh2.ne.cnh1) then                              ! no further particles in cell

               spx_all=0
               do sp1=0,spx_nsp
                  spx_all=spx_all+spx_np(sp1)
               enddo

               if (spx_all.gt.1) then                          ! more than one particle in cell

               
c PERFORM INTRA-SPECIES COLLISIONS

c                  goto 1111                                     ! use to turn same particle collisions off
                  do sp1=0,spx_nsp

                     if (spx_np(sp1).gt.1) then

                        rho=0.0
                        do j=0,spx_np(sp1)-1                      ! over all particles of sort sp1
                           wni=p_niloc(11*spx_list(sp1,j)+10)     ! weights of particles of sort sp1 in cell
                           rho=rho+wni                            ! cumulative weight of all particles of sort sp1 in cell
                        enddo
                        rho=cori*n0*rho                           ! density of specie sp1                

c                        write(6,*) 'rho:',rho,spx_np(sp1)
c                        rho=cori*n0*spx_np(sp1)
c                        write(6,*) 'rho:',rho,spx_np(sp1)

                        if (mod(spx_np(sp1),2).eq.1) then         ! odd number of particles

                           call bc(rho,lnL,spx_list(sp1,0), 
     &                                     spx_list(sp1,1),1)

                           coll_count(sp1,sp1)=
     &                        coll_count(sp1,sp1)+coll_add                   ! coll counter
                           if (nudt.gt.nudt_max(sp1,sp1)) 
     &                        nudt_max(sp1,sp1)=nudt                         ! coll freq min
                           if (nudt.lt.nudt_min(sp1,sp1)) 
     &                        nudt_min(sp1,sp1)=nudt                         ! coll freq max
                           nudt_mean(sp1,sp1)=
     &                        nudt_mean(sp1,sp1)+nudt                        ! coll freq mean

                           call bc(rho,lnL,spx_list(sp1,0), 
     &                                     spx_list(sp1,2),1)

                           coll_count(sp1,sp1)=
     &                        coll_count(sp1,sp1)+coll_add                   ! coll counter
                           if (nudt.gt.nudt_max(sp1,sp1)) 
     &                        nudt_max(sp1,sp1)=nudt                         ! coll freq min
                           if (nudt.lt.nudt_min(sp1,sp1)) 
     &                        nudt_min(sp1,sp1)=nudt                         ! coll freq max
                           nudt_mean(sp1,sp1)=
     &                        nudt_mean(sp1,sp1)+nudt                        ! coll freq mean

                           call bc(rho,lnL,spx_list(sp1,1), 
     &                                     spx_list(sp1,2),1)

                           coll_count(sp1,sp1)=
     &                        coll_count(sp1,sp1)+coll_add                   ! coll counter
                           if (nudt.gt.nudt_max(sp1,sp1)) 
     &                        nudt_max(sp1,sp1)=nudt                         ! coll freq min
                           if (nudt.lt.nudt_min(sp1,sp1)) 
     &                        nudt_min(sp1,sp1)=nudt                         ! coll freq max
                           nudt_mean(sp1,sp1)=
     &                        nudt_mean(sp1,sp1)+nudt                        ! coll freq mean

                           if (spx_np(sp1).gt.3) then
                              do i=3,spx_np(sp1)-1,2
                                 call bc(rho,lnL,spx_list(sp1,i), 
     &                                           spx_list(sp1,i+1),0)

                                 coll_count(sp1,sp1)=
     &                              coll_count(sp1,sp1)+coll_add             ! coll counter
                                 if (nudt.gt.nudt_max(sp1,sp1)) 
     &                              nudt_max(sp1,sp1)=nudt                   ! coll freq min
                                 if (nudt.lt.nudt_min(sp1,sp1)) 
     &                              nudt_min(sp1,sp1)=nudt                   ! coll freq max
                                 nudt_mean(sp1,sp1)=
     &                              nudt_mean(sp1,sp1)+nudt                  ! coll freq mean

                              enddo
                           endif

                        else                      ! even number of particles

                           do i=0,spx_np(sp1)-1,2
                              call bc(rho,lnL,spx_list(sp1,i),
     &                                        spx_list(sp1,i+1),0)

                              coll_count(sp1,sp1)=
     &                           coll_count(sp1,sp1)+coll_add                ! coll counter
                              if (nudt.gt.nudt_max(sp1,sp1)) 
     &                           nudt_max(sp1,sp1)=nudt                      ! coll freq min
                              if (nudt.lt.nudt_min(sp1,sp1)) 
     &                           nudt_min(sp1,sp1)=nudt                      ! coll freq max
                              nudt_mean(sp1,sp1)=
     &                           nudt_mean(sp1,sp1)+nudt                     ! coll freq mean

                           enddo

                        endif
                     endif 
                  enddo
c 1111             continue

c PERFORM INTER-SPECIES COLLISIONS 

                  
c                  goto 2222                                     ! use to turn inter-particle collisions off
                  do sp1=0,spx_nsp
                  do sp2=sp1+1,spx_nsp
 
                     if ((spx_np(sp1).gt.0).and.(spx_np(sp2).gt.0)) then

c                        rho=cori*n0*min(spx_np(sp1),spx_np(sp2))     ! min density of species in the pair
c                        write(6,*) 'rho1:',rho,spx_np(sp1),spx_np(sp2)

                        if (spx_np(sp1).eq.spx_np(sp2)) then             ! N(SP1)==N(SP2)

                           rho=0.0
                           do j=0,min(spx_np(sp1),spx_np(sp2))-1     ! over all particles of sort sp1
                              wni=p_niloc(11*spx_list(sp2,j)+10)     ! weights of particles of sort sp1 in cell
                              rho=rho+wni                            ! cumulative weight of all particles of sort sp1 in cell
                           enddo
                           rho=cori*n0*rho                           ! min density of specie sp1 and sp2                
c                           write(6,*) 'rho2:',rho

                           do i=0,spx_np(sp1)-1
                              call bc(rho,lnL,spx_list(sp1,i),
     &                                        spx_list(sp2,i),0)

                              coll_count(sp1,sp2)=
     &                           coll_count(sp1,sp2)+coll_add                 ! coll counter
                              if (nudt.gt.nudt_max(sp1,sp2)) 
     &                           nudt_max(sp1,sp2)=nudt                       ! coll freq min
                              if (nudt.lt.nudt_min(sp1,sp2)) 
     &                           nudt_min(sp1,sp2)=nudt                       ! coll freq max
                              nudt_mean(sp1,sp2)=
     &                           nudt_mean(sp1,sp2)+nudt                      ! coll freq mean

                           enddo

                        else 

                           if (spx_np(sp1).gt.spx_np(sp2)) then          ! N(sp1)>N(sp2)

                              rho=0.0
                              do j=0,min(spx_np(sp1),spx_np(sp2))-1     ! over all particles of sort sp1
                                 wni=p_niloc(11*spx_list(sp2,j)+10)     ! weights of particles of sort sp1 in cell
                                 rho=rho+wni                            ! cumulative weight of all particles of sort sp1 in cell
                              enddo
                              rho=cori*n0*rho                           ! min density of specie sp1 and sp2                
c                              write(6,*) 'rho3:',rho

                              intg=int(1.0d0*spx_np(sp1)/spx_np(sp2))
                              rest=1.0d0*spx_np(sp1)/spx_np(sp2)-intg

                              limit=int(rest*spx_np(sp2)+0.5d0) 
                              do i=0,limit-1
                                 do j=i*(intg+1),(i+1)*(intg+1)-1
                                    call bc(rho,lnL,spx_list(sp1,j), 
     &                                              spx_list(sp2,i),0)

                                    coll_count(sp1,sp2)=
     &                                 coll_count(sp1,sp2)+coll_add           ! coll counter
                                    if (nudt.gt.nudt_max(sp1,sp2)) 
     &                                 nudt_max(sp1,sp2)=nudt                 ! coll freq min
                                    if (nudt.lt.nudt_min(sp1,sp2)) 
     &                                 nudt_min(sp1,sp2)=nudt                 ! coll freq max
                                    nudt_mean(sp1,sp2)=
     &                                 nudt_mean(sp1,sp2)+nudt                ! coll freq mean

                                 enddo
                              enddo

                              j=limit*(intg+1)
                              do i=limit,spx_np(sp2)-1
                                 do k=0,intg-1
                                    call bc(rho,lnL,spx_list(sp1,j), 
     &                                              spx_list(sp2,i),0)

                                    coll_count(sp1,sp2)=
     &                                 coll_count(sp1,sp2)+coll_add           ! coll counter
                                    if (nudt.gt.nudt_max(sp1,sp2)) 
     &                                 nudt_max(sp1,sp2)=nudt                 ! coll freq min
                                    if (nudt.lt.nudt_min(sp1,sp2)) 
     &                                 nudt_min(sp1,sp2)=nudt                 ! coll freq max
                                    nudt_mean(sp1,sp2)=
     &                                 nudt_mean(sp1,sp2)+nudt                ! coll freq mean

                                    j=j+1
                                 enddo
                              enddo
                               
                           else                                          ! N(sp1)<N(sp2)

                              rho=0.0
                              do j=0,min(spx_np(sp1),spx_np(sp2))-1     ! over all particles of sort sp1
                                 wni=p_niloc(11*spx_list(sp1,j)+10)     ! weights of particles of sort sp1 in cell
                                 rho=rho+wni                            ! cumulative weight of all particles of sort sp1 in cell
                              enddo
                              rho=cori*n0*rho                           ! min density of specie sp1 and sp2                
c                              write(6,*) 'rho4:',rho

                              intg=int(1.0d0*spx_np(sp2)/spx_np(sp1))
                              rest=1.0d0*spx_np(sp2)/spx_np(sp1)-intg

                              limit=int(rest*spx_np(sp1)+0.5d0) 
                              do i=0,limit-1
                                 do j=i*(intg+1),(i+1)*(intg+1)-1
                                    call bc(rho,lnL,spx_list(sp2,j), 
     &                                              spx_list(sp1,i),0)

                                    coll_count(sp1,sp2)=
     &                                 coll_count(sp1,sp2)+coll_add           ! coll counter
                                    if (nudt.gt.nudt_max(sp1,sp2)) 
     &                                 nudt_max(sp1,sp2)=nudt                 ! coll freq min
                                    if (nudt.lt.nudt_min(sp1,sp2)) 
     &                                 nudt_min(sp1,sp2)=nudt                 ! coll freq max
                                    nudt_mean(sp1,sp2)=
     &                                 nudt_mean(sp1,sp2)+nudt                ! coll freq mean

                                 enddo
                              enddo

                              j=limit*(intg+1)
                              do i=limit,spx_np(sp1)-1
                                 do k=0,intg-1
                                    call bc(rho,lnL,spx_list(sp2,j), 
     &                                              spx_list(sp1,i),0)

                                    coll_count(sp1,sp2)=
     &                                 coll_count(sp1,sp2)+coll_add           ! coll counter
                                    if (nudt.gt.nudt_max(sp1,sp2)) 
     &                                 nudt_max(sp1,sp2)=nudt                 ! coll freq min
                                    if (nudt.lt.nudt_min(sp1,sp2)) 
     &                                 nudt_min(sp1,sp2)=nudt                 ! coll freq max
                                    nudt_mean(sp1,sp2)=
     &                                 nudt_mean(sp1,sp2)+nudt                ! coll freq mean

                                    j=j+1
                                 enddo
                              enddo
                               
                           endif
                        endif
                     endif

                  enddo
                  enddo
c 2222             continue


c END COLLISIONS

          
               endif
             
               spx_list=0       ! reset list
               spx_np=0         ! reset counters
               
            endif
         enddo
         

c OUTPUT OF CRITICAL COLLISION PARAMETERS


         if (n.eq.nprc) then
            nprc=nprc+dnprc

            call SERV_labelgen(mpe,node)
            call SERV_labelgen(n,label)

            open(11,file=trim(data_out)
     &           //'/'//node//'collision'//label,
     &           access='sequential',form='formatted')

            write(11,*) 'COLLISIONS PARAMETERS at timestep:',n
            do sp1=0,spx_nsp
               do sp2=sp1,spx_nsp
                  coll_tot=coll_tot+coll_count(sp1,sp2)
                  write(11,*) ' '
                  write(11,*) 'particle pair:',sp1,sp2
                  write(11,*) 'coll_count:',coll_count(sp1,sp2)
                  write(11,*) 'nudt_min:',nudt_min(sp1,sp2)
                  write(11,*) 'nudt_max:',nudt_max(sp1,sp2)
                  write(11,*) 'nudt_mean:',nudt_mean(sp1,sp2)
     &                         /(coll_count(sp1,sp2)+1.0d-10)
                  write(11,*) ' '
               enddo
            enddo
            write(11,*) 'total number of collisions:',coll_tot

            close(11)      

         endif
      endif


      deallocate(spx_list)
      deallocate(spx_np)

      deallocate(coll_count)
      deallocate(nudt_min)
      deallocate(nudt_max)
      deallocate(nudt_mean)


      end subroutine PIC_bin_coll






C PROCESS THE BINARY COLLISION BETWEEN ANY TWO GIVEN PARTICLES


      subroutine bc(p_niloc, rho,lnL,indl,indm,nu_half)

      use PIC_bin_collision
      use VLA_variables, only: pi
      use PSC_domain, only: npe, mpe
      use PSC_params
      use PSC_patch

      implicit none

      real(kind=8) :: p_niloc(0:*)

      integer :: indl,indm,nu_half

      real(kind=8) :: pn1,pn2,pn3,pn4
      real(kind=8) :: p01,p02,p03,p04,pc01,pc02,pc03,pc04
      real(kind=8) :: px1,py1,pz1,pcx1,pcy1,pcz1
      real(kind=8) :: px2,py2,pz2,pcx2,pcy2,pcz2
      real(kind=8) :: px3,py3,pz3,pcx3,pcy3,pcz3
      real(kind=8) :: px4,py4,pz4,pcx4,pcy4,pcz4
      real(kind=8) :: h1,h2,h3,h4,ppc,qqc,ss
      real(kind=8) :: rho,lnL

      real(kind=8) :: m1,m2,m3,m4
      real(kind=8) :: q1,q2,q3,q4
      real(kind=8) :: w1,w2,w3,w4,ww
      real(kind=8) :: vcx,vcy,vcz
      real(kind=8) :: bet,gam

      real(kind=8) :: psi,nu
      real(kind=8) :: nx,ny,nz,nnorm
      real(kind=8) :: nn1,nx1,ny1,nz1
      real(kind=8) :: nn2,nx2,ny2,nz2
      real(kind=8) :: nn3,nx3,ny3,nz3
      real(kind=8) :: vcx1,vcy1,vcz1
      real(kind=8) :: vcx2,vcy2,vcz2
      real(kind=8) :: vcn,vcxr,vcyr,vczr,vcr
      real(kind=8) :: m12,q12
      real(kind=8) :: ran1,ran2,ran3
      real(kind=8) :: pp1,pp2

      real(kind=8),allocatable,dimension(:) :: rndmv


      allocate(rndmv(1:3*npe))

      px1=p_niloc(11*indl+3)
      py1=p_niloc(11*indl+4)
      pz1=p_niloc(11*indl+5)
      q1=p_niloc(11*indl+6)
      m1=p_niloc(11*indl+7)
      w1=p_niloc(11*indl+10)
      px2=p_niloc(11*indm+3)
      py2=p_niloc(11*indm+4)
      pz2=p_niloc(11*indm+5)
      q2=p_niloc(11*indm+6)
      m2=p_niloc(11*indm+7)
      w2=p_niloc(11*indm+10)


      px1=m1*px1
      py1=m1*py1
      pz1=m1*pz1
      px2=m2*px2
      py2=m2*py2
      pz2=m2*pz2


c determine absolute value of pre-collision momentum in cm-frame


      p01=dsqrt(m1*m1+px1*px1+py1*py1+pz1*pz1)
      p02=dsqrt(m2*m2+px2*px2+py2*py2+pz2*pz2)
      h1=p01*p02-px1*px2-py1*py2-pz1*pz2
      ss=m1*m1+m2*m2+2.0*h1
      h2=ss-m1*m1-m2*m2
      h3=(h2*h2-4.0*m1*m1*m2*m2)/(4.0*ss)
      ppc=dsqrt(abs(h3))


c determine cm-velocity


      vcx=(px1+px2)/(p01+p02)
      vcy=(py1+py2)/(p01+p02)
      vcz=(pz1+pz2)/(p01+p02)
      
      nnorm=dsqrt(vcx*vcx+vcy*vcy+vcz*vcz)
      if (nnorm>0.0d0) then
         nx=vcx/nnorm
         ny=vcy/nnorm
         nz=vcz/nnorm
      else
         nx=0.0d0
         ny=0.0d0
         nz=0.0d0
      endif
      bet=nnorm
      gam=1.0/dsqrt(1.0-bet*bet)


c determine pre-collision momenta in cm-frame
      
      
      pn1=px1*nx+py1*ny+pz1*nz
      pn2=px2*nx+py2*ny+pz2*nz
      pc01=dsqrt(m1*m1+ppc*ppc)
      pcx1=px1+(gam-1.0)*pn1*nx-gam*vcx*p01
      pcy1=py1+(gam-1.0)*pn1*ny-gam*vcy*p01
      pcz1=pz1+(gam-1.0)*pn1*nz-gam*vcz*p01
      pc02=dsqrt(m2*m2+ppc*ppc)
      pcx2=px2+(gam-1.0)*pn2*nx-gam*vcx*p02
      pcy2=py2+(gam-1.0)*pn2*ny-gam*vcy*p02
      pcz2=pz2+(gam-1.0)*pn2*nz-gam*vcz*p02


c introduce right-handed coordinate system
      
      
      nn1=dsqrt(pcx1*pcx1+pcy1*pcy1+pcz1*pcz1)
      nn2=dsqrt(pcx1*pcx1+pcy1*pcy1)
      nn3=nn1*nn2

      if (nn2.ne.0.0d0) then
         nx1=pcx1/nn1
         ny1=pcy1/nn1
         nz1=pcz1/nn1
         
         nx2=pcy1/nn2
         ny2=-pcx1/nn2
         nz2=0.0
         
         nx3=-pcx1*pcz1/nn3
         ny3=-pcy1*pcz1/nn3
         nz3=nn2*nn2/nn3
      else
         nx1=0.0d0
         ny1=0.0d0
         nz1=1.0d0

         nx2=0.0d0
         ny2=1.0d0
         nz2=0.0d0

         nx3=1.0d0
         ny3=0.0d0
         nz3=0.0d0
      endif


c determine relative particle velocity in cm-frame


      vcx1=pcx1/pc01
      vcy1=pcy1/pc01
      vcz1=pcz1/pc01
      vcx2=pcx2/pc02
      vcy2=pcy2/pc02
      vcz2=pcz2/pc02
      
      vcn=1.0d0/(1.0d0-(vcx1*vcx2+vcy1*vcy2+vcz1*vcz2))
      vcxr=vcn*(vcx1-vcx2)
      vcyr=vcn*(vcy1-vcy2)
      vczr=vcn*(vcz1-vcz2)
      vcr=max(1.0d-20,dsqrt(vcxr*vcxr+vcyr*vcyr+vczr*vczr))

      m3=m1
      m4=m2
      q3=q1
      q4=q2
      w3=w1
      w4=w2


c determine absolute value of post-collision momentum in cm-frame


      if (ss.ge.(m3+m4)*(m3+m4)) then

         h2=ss-m3*m3-m4*m4
         h3=(h2*h2-4.0*m3*m3*m4*m4)/(4.0*ss)
         qqc=dsqrt(abs(h3))
         m12=m1*m2/(m1+m2)
         q12=q1*q2

         if (q12*q12.gt.0.0) then

c Spitzer collision frequency
c v_{ei}=q1^2*q2^2*ni*lnL/(4pi*eps0^2*mr^2*vr^3)
c v_{ei}=3.3*10^7*[(q1/e)^2*(q2/e)^2*(ni/10^27)*lnL]/[(mr/me)^2*(vr/c)^3]
c 2.9849366920d-20=qq**4/(4.0*pi*eps0**2*mm**2*c**3)

            coll_add=1.0d0
            nudt=2.9849366920d-20
            nudt=nudt*lnL*q12*q12*rho*dt/(wl*m12*m12*vcr*vcr*vcr)
            if (nu_half.eq.0) nudt=1.0d0*nudt
            if (nu_half.eq.1) nudt=0.5d0*nudt

c event generator of angles for post collision vectors

            call random_number(rndmv)
            ran1=max(1.0d-20,rndmv(3*mpe+1))
            ran2=max(1.0d-20,rndmv(3*mpe+2))
            ran3=max(1.0d-20,rndmv(3*mpe+3))

            ww=max(w1,w2)
            if (ran1.le.w2/ww) then
               pp1=1.0
            else
               pp1=0.0
            endif
            if (ran2.le.w1/ww) then
               pp2=1.0
            else
               pp2=0.0
            endif

            nu=2.0d0*pi*ran2

            if(nudt<0.5d0) then                                     ! small angle collision
               psi=2.0*atan(sqrt(-0.5*nudt*log(1.0-ran3)))
            else
               psi=acos(1.0-2.0*ran3)                               ! isotropic angles
            endif

         endif


c         if (q12*q12.eq.0.0) then

c Hard sphere collision frequency
c sigma_0=1.0d-20 m^2/4*pi

c            coll_add=1.0d0
c            nudt=1.0d-20
c            nudt=nudt*rho*vcr*cc*dt/wl

c event generator of angles for post collision vectors

c            call random_number(rndmv)
c            ran1=max(1.0d-20,rndmv(3*mpe+1))
c            ran2=max(1.0d-20,rndmv(3*mpe+2))
c            ran3=max(1.0d-20,rndmv(3*mpe+3))

c            ww=max(w1,w2)
c            if (ran1.le.w2/ww) then
c               pp1=1.0
c            else
c               pp1=0.0
c            endif
c            if (ran2.le.w1/ww) then
c               pp2=1.0
c            else
c               pp2=0.0
c            endif

c            nu=2.0d0*pi*ran2
c            psi=acos(1.0-2.0*ran3)                                 ! isotropic angle distribution for neutrals

c         endif

      else
         coll_add=0.0d0
         nudt=0.0d0
         return
      endif
      

c determine post-collision momentum in cm-frame


      h1=cos(psi)
      h2=sin(psi)
      h3=sin(nu)                        
      h4=cos(nu)

      if (pp1.eq.1.0) then
         pc03=dsqrt(m3*m3+qqc*qqc)
         pcx3=qqc*(h1*nx1+h2*h3*nx2+h2*h4*nx3)
         pcy3=qqc*(h1*ny1+h2*h3*ny2+h2*h4*ny3)
         pcz3=qqc*(h1*nz1+h2*h3*nz2+h2*h4*nz3)
      else
         pc03=pc01
         pcx3=pcx1
         pcy3=pcy1
         pcz3=pcz1
      endif
      if (pp2.eq.1.0) then
         pc04=dsqrt(m4*m4+qqc*qqc)
         pcx4=-qqc*(h1*nx1+h2*h3*nx2+h2*h4*nx3)
         pcy4=-qqc*(h1*ny1+h2*h3*ny2+h2*h4*ny3)
         pcz4=-qqc*(h1*nz1+h2*h3*nz2+h2*h4*nz3)
      else
         pc04=pc02
         pcx4=pcx2
         pcy4=pcy2
         pcz4=pcz2
      endif      


c determine post-collision momentum in lab-frame

      
      pn3=pcx3*nx+pcy3*ny+pcz3*nz
      pn4=pcx4*nx+pcy4*ny+pcz4*nz
      p03=gam*(pc03+bet*pn3)
      px3=pcx3+(gam-1.0)*pn3*nx+gam*vcx*pc03
      py3=pcy3+(gam-1.0)*pn3*ny+gam*vcy*pc03
      pz3=pcz3+(gam-1.0)*pn3*nz+gam*vcz*pc03
      p04=gam*(pc04+bet*pn4)
      px4=pcx4+(gam-1.0)*pn4*nx+gam*vcx*pc04
      py4=pcy4+(gam-1.0)*pn4*ny+gam*vcy*pc04
      pz4=pcz4+(gam-1.0)*pn4*nz+gam*vcz*pc04
     
      px3=px3/m3
      py3=py3/m3
      pz3=pz3/m3
      px4=px4/m4
      py4=py4/m4
      pz4=pz4/m4

      p_niloc(11*indl+3)=px3
      p_niloc(11*indl+4)=py3
      p_niloc(11*indl+5)=pz3
      p_niloc(11*indl+7)=m3
      p_niloc(11*indl+10)=w3
      p_niloc(11*indm+3)=px4
      p_niloc(11*indm+4)=py4
      p_niloc(11*indm+5)=pz4
      p_niloc(11*indm+7)=m4
      p_niloc(11*indm+10)=w4


      deallocate(rndmv)


      end subroutine bc
