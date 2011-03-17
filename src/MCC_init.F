c     INITIALIZE IMPACT IONIZATION PACKAGE BEFORE CALLING MCC_impact:

      subroutine INIT_MCC

      use MCC_variables
      use PSC_patch

      implicit none

      integer :: mat
      integer :: cs
      integer :: NDAT
      
      real(kind=8) :: ek,ekmax,sv,va


      NMAT          =1                                  ! number of materials
                                                        ! in case of NMAT>1, change table allocation
      NCS           =2                                  ! number of allowed charge states>0

      NPMAX         =3*(NCS+1)*nicell                   ! conservative estimate for max number of particles in cell

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
      NDAT          =125                                ! max number of impact x-section data points
      me            =5.11d5                             ! electron rest mass in eV

      allocate(matname(1:NMAT))



      allocate(n_xstable(1:NMAT)) 
      allocate(xstable_n(1:NMAT, 0:NCS-1))
      allocate(xstable_t(1:NMAT, 0:NCS-1))
      allocate(max_sigmav(1:NMAT, 0:NCS-1))
      allocate(xstable_e(1:NMAT, 0:NCS-1, 0:NDAT))
      allocate(xstable_s(1:NMAT, 0:NCS-1, 0:NDAT))
       

c     INITIALIZE IMPACT TABLES

      n_xstable =(/ 2, 2 /)                                    ! number of tables == charge states>0
      matname(:)=(/'HELIUM','CARBON'/)                         ! names of materials

c     DEFAULT VALUES FOR ELECTRONS

      mat=0
      mpart(mat)=1.0                                           ! electron mass

c     IMPACT IONIZATION OF HELIUM

      mat=1
      if(mat.gt.NMAT) then
         write(*,*) 'material not allowed'
         stop
      endif

      mpart(mat)         =7344.0                               ! use ion mass to recognize material
      xstable_n(mat,:)   =(/ 122, 16 /)                        ! number of array elements < NDAT
      write(6,*) 'xstable_n:',xstable_n(1,0),xstable_n(1,1)
      xstable_t(mat,:)   =(/ 24.59d0 , 55.0d0  /)              ! E_threshold in eV  
                                                               ! sigma(Ekin): [m^2] ([eV])
      
      cs=0                                               ! He, Z=0-->1
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      xstable_e(mat,cs,:)=(/ 24.59, 26., 26.5, 26.6, 27., 27.5, 
     c     27.6, 28., 28.5, 28.6, 29., 29.5, 29.6, 30., 30.5, 30.6, 31., 
     c     31.5, 32., 32.1, 32.5, 33., 33.5, 33.6, 34., 34.5, 35., 35.5, 
     c     36., 37., 38., 38.6, 39., 40., 41., 42., 43., 
     c     43.6, 44., 45., 46., 47., 48., 48.6, 49., 50., 51.,  
     c     52., 53.6, 54., 56., 58., 58.6, 60., 65., 67.,
     c     68.6, 70., 75., 78.6, 80., 85., 88.6, 90., 90.2, 95., 95.2, 
     c     100., 
     c     105., 110., 115., 118., 120., 130., 135., 140., 145., 147., 
     c     150., 160., 
     c     170., 180., 190., 200., 220., 225., 250., 275., 280., 300., 
     c     325., 350., 
     c     375., 400., 430., 450., 500., 550., 570., 600., 650., 700., 
     c     750., 870., 
     c     1000., 1150., 1320., 1520., 1750., 2010., 2300., 2650., 
     c     3000., 3500., 
     c     4000., 4600., 5300., 6100., 7000., 8000., 9000., 10000. /)
      xstable_s(mat,cs,:)=(/0., 1.9e-22, 2.6e-22, 2.8e-22, 3.3e-22,  
     c     4e-22, 4.1e-22, 4.7e-22, 5.4e-22, 5.5e-22, 6e-22, 6.7e-22, 
     c     6.8e-22, 
     c     7.3e-22, 8e-22, 8.1e-22, 8.6e-22, 9.2e-22, 9.8e-22, 1e-21, 
     c     1.04e-21, 1.1e-21, 1.16e-21, 1.17e-21, 1.22e-21, 1.27e-21, 
     c     1.33e-21, 1.38e-21, 1.43e-21, 1.54e-21, 1.63e-21, 1.69e-21, 
     c     1.73e-21, 1.82e-21, 1.9e-21, 1.98e-21, 2.06e-21, 2.11e-21, 
     c     2.14e-21, 2.21e-21, 2.28e-21, 2.34e-21, 2.41e-21, 2.44e-21, 
     c     2.47e-21, 2.52e-21, 2.58e-21, 2.63e-21, 2.71e-21, 2.73e-21, 
     c     2.82e-21, 2.9e-21, 2.92e-21, 2.98e-21, 3.14e-21, 3.19e-21, 
     c     3.23e-21, 3.26e-21, 3.36e-21, 3.42e-21, 3.44e-21, 3.5e-21, 
     c     3.53e-21, 3.54e-21, 3.55e-21, 3.58e-21, 3.58e-21, 3.6e-21, 
     c     3.62e-21, 3.62e-21, 3.62e-21, 3.62e-21, 3.62e-21, 3.6e-21, 
     c     3.59e-21, 3.57e-21, 3.55e-21, 3.54e-21, 3.53e-21, 3.49e-21, 
     c     3.44e-21, 3.39e-21, 3.33e-21, 3.28e-21, 3.17e-21, 3.15e-21, 
     c     3.02e-21, 2.9e-21, 2.87e-21, 2.78e-21, 2.67e-21, 2.57e-21, 
     c     2.48e-21, 2.39e-21, 2.3e-21, 2.23e-21, 2.1e-21, 1.97e-21, 
     c     1.93e-21, 1.87e-21, 1.77e-21, 1.68e-21, 1.61e-21, 1.45e-21, 
     c     1.31e-21, 1.18e-21, 1.06e-21, 9.5e-22, 8.5e-22, 7.6e-22, 
     c     6.8e-22, 6.1e-22, 5.5e-22, 4.8e-22, 4.3e-22, 3.8e-22, 
     c     3.4e-22, 3e-22, 2.7e-22, 2.4e-22, 2.2e-22, 2e-22 /)

      cs=1                                                              ! He, Z=1-->2
      if(cs.gt.NCS-1) then 
         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
         stop
      endif
      xstable_e(mat,cs,:)=(/ 55., 70., 80., 90., 100., 150., 200., 300., 
     c     400., 500., 
     c     700., 1000., 2000., 3000., 10000., 100000. /)
      xstable_s(mat,cs,:)=(/ 0., 1.5e-22, 2.75e-22, 3.25e-22, 3.55e-22, 
     c     4.5e-22, 4.5e-22, 4.25e-22, 3.9e-22, 3.5e-22, 3e-22, 
     c     2.15e-22, 1.25e-22, 8e-23, 4e-23, 5e-24 /)





c     IMPACT IONIZATION OF CARBON

c      mat=2
c      if(mat.gt.NMAT) then
c         write(*,*) 'material not allowed'
c         stop
c      endif

c      mpart(mat)  =27000.0 

c      n_xstable(mat)   =2                                    ! number of tables == charge states>0
c      xstable_n(mat,:) =(/ 122, 16 /)                        ! number of array elements < NDAT
c      xstable_t(mat,:) =(/ 24.59d0 , 55.0d0  /)              ! E_threshold in eV  
                                                             ! sigma(Ekin): [m^2] ([eV])

c      cs=0                                                   ! He, Z=0-->1
c      if(cs.gt.NCS-1) then 
c         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
c         stop
c      endif
c      xstable_e(mat,cs,:)=(/ 24.59, 26., 26.5, 26.6, 27., 27.5, 
c     c     27.6, 28., 28.5, 28.6, 29., 29.5, 29.6, 30., 30.5, 30.6, 31., 
c     c     31.5, 32., 32.1, 32.5, 33., 33.5, 33.6, 34., 34.5, 35., 35.5, 
c     c     36., 37., 38., 38.6, 39., 40., 41., 42., 43., 
c     c     43.6, 44., 45., 46., 47., 48., 48.6, 49., 50., 51.,  
c     c     52., 53.6, 54., 56., 58., 58.6, 60., 65., 67.,
c     c     68.6, 70., 75., 78.6, 80., 85., 88.6, 90., 90.2, 95., 95.2, 
c     c     100., 
c     c     105., 110., 115., 118., 120., 130., 135., 140., 145., 147., 
c     c     150., 160., 
c     c     170., 180., 190., 200., 220., 225., 250., 275., 280., 300., 
c     c     325., 350., 
c     c     375., 400., 430., 450., 500., 550., 570., 600., 650., 700., 
c     c     750., 870., 
c     c     1000., 1150., 1320., 1520., 1750., 2010., 2300., 2650., 
c     c     3000., 3500., 
c     c     4000., 4600., 5300., 6100., 7000., 8000., 9000., 10000. /)
c      xstable_s(mat,cs,:)=(/0., 1.9e-22, 2.6e-22, 2.8e-22, 3.3e-22, 
c     c     4e-22, 
c     c     4.1e-22, 4.7e-22, 5.4e-22, 5.5e-22, 6e-22, 6.7e-22, 
c     c     6.8e-22, 
c     c     7.3e-22, 8e-22, 8.1e-22, 8.6e-22, 9.2e-22, 9.8e-22, 1e-21, 
c     c     1.04e-21, 1.1e-21, 1.16e-21, 1.17e-21, 1.22e-21, 1.27e-21, 
c     c     1.33e-21, 1.38e-21, 1.43e-21, 1.54e-21, 1.63e-21, 1.69e-21, 
c     c     1.73e-21, 1.82e-21, 1.9e-21, 1.98e-21, 2.06e-21, 2.11e-21, 
c     c     2.14e-21, 2.21e-21, 2.28e-21, 2.34e-21, 2.41e-21, 2.44e-21, 
c     c     2.47e-21, 2.52e-21, 2.58e-21, 2.63e-21, 2.71e-21, 2.73e-21, 
c     c     2.82e-21, 2.9e-21, 2.92e-21, 2.98e-21, 3.14e-21, 3.19e-21, 
c     c     3.23e-21, 3.26e-21, 3.36e-21, 3.42e-21, 3.44e-21, 3.5e-21, 
c     c     3.53e-21, 3.54e-21, 3.55e-21, 3.58e-21, 3.58e-21, 3.6e-21, 
c     c     3.62e-21, 3.62e-21, 3.62e-21, 3.62e-21, 3.62e-21, 3.6e-21, 
c     c     3.59e-21, 3.57e-21, 3.55e-21, 3.54e-21, 3.53e-21, 3.49e-21, 
c     c     3.44e-21, 3.39e-21, 3.33e-21, 3.28e-21, 3.17e-21, 3.15e-21, 
c     c     3.02e-21, 2.9e-21, 2.87e-21, 2.78e-21, 2.67e-21, 2.57e-21, 
c     c     2.48e-21, 2.39e-21, 2.3e-21, 2.23e-21, 2.1e-21, 1.97e-21, 
c     c     1.93e-21, 1.87e-21, 1.77e-21, 1.68e-21, 1.61e-21, 1.45e-21, 
c     c     1.31e-21, 1.18e-21, 1.06e-21, 9.5e-22, 8.5e-22, 7.6e-22, 
c     c     6.8e-22, 6.1e-22, 5.5e-22, 4.8e-22, 4.3e-22, 3.8e-22, 
c     c     3.4e-22, 3e-22, 2.7e-22, 2.4e-22, 2.2e-22, 2e-22 /)

c      cs=1                                                              ! He, Z=1-->2
c      if(cs.gt.NCS-1) then 
c         write(*,*) 'ERROR IN INITIALZING IMPACT VARIABLES'
c         stop
c      endif
c      xstable_e(mat,cs,:)=(/ 55., 70., 80., 90., 100., 150., 200., 300., 
c     c     400., 500., 
c     c     700., 1000., 2000., 3000., 10000., 100000. /)
c      xstable_s(mat,cs,:)=(/ 0., 1.5e-22, 2.75e-22, 3.25e-22, 3.55e-22, 
c     c     4.5e-22, 4.5e-22, 4.25e-22, 3.9e-22, 3.5e-22, 3e-22, 
c     c     2.15e-22, 1.25e-22, 8e-23, 4e-23, 5e-24 /)

      

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
