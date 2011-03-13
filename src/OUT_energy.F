c THIS SUBROUTINE WRITES DETAILED KINETIC TEMPERATURES FOR ELECTRONS AND IONS
c by A.Kemp 11/15/2005

      subroutine OUT_energy(niloc, p_niloc)

      use VLA_variables
      use OUT_fields
      use PSC_params
      use PSC_patch

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)

      character*(5) node, label

      character*100 :: fname='energy.data'
      character*100 :: fnames='-spectrum.data'

      real(kind=8) :: te,ti,ekini,ekine,ekin
      real(kind=8) :: ekinl,ekinx,dekin
      real(kind=8) :: ikinl,ikinx,dikin
      real(kind=8) :: xne,xni,qni,mni,wni,ux,uy,uz
      integer      :: l,bek,ntot

      real(kind=8),allocatable,dimension(:) :: enwe,enwi


      ekine=0.0
      ekini=0.0
      xne=0.0
      xni=0.0

      ekinl=0.0
      ekinx=3.0 ! [MeV]
      dekin=ekinx-ekinl

      ikinl=0.000000
      ikinx=0.020000 ! MeV
      dikin=ikinx-ikinl

      allocate(enwe(0:399))
      allocate(enwi(0:399))
      enwe=0                    ! reset array
      ntot=0

      do l=1,niloc
	
	qni=p_niloc(11*l+6)
	mni=p_niloc(11*l+7)

        ux=p_niloc(11*l+3)
	uy=p_niloc(11*l+4)
	uz=p_niloc(11*l+5)

 	wni=p_niloc(11*l+10)

	ekin=0.511733d0*wni*mni*(dsqrt(1.d0+ux*ux+uy*uy+uz*uz)-1.d0)   ! in MeV

	if (qni==-1.d0) then 

	   ekine=ekine+ekin
           xne=xne+1.d0*wni

           if( ekin>=ekinl.AND.ekin<=ekinx) then
              bek= int((ekin-ekinl)*399.d0/dekin+0.5)
              if(bek>=0.AND.bek<=399) then 
                 enwe(bek)=enwe(bek)+1 
                 ntot=ntot+1
              else 
                 write(*,*) 'bek out of range'
                 call MPI_FINALIZE(info)
                 stop
              endif
           else
              write(*,*) 'ekin out of range', ekin
              call MPI_FINALIZE(info)
              stop
           endif

        else 

 	   ekini=ekini+ekin
	   xni=xni+1.d0*wni

           if( ekin>=ikinl.AND.ekin<=ikinx) then
              bek= int((ekin-ikinl)*399.d0/dikin+0.5)
              if(bek>=0.AND.bek<=399) then 
                 enwi(bek)=enwi(bek)+1 
              else 
                 write(*,*) 'bek-i out of range'
                 call MPI_FINALIZE(info)
                 stop
              endif
           else
              write(*,*) 'ekin-i out of range', ekin
              call MPI_FINALIZE(info)
              stop
           endif
 	endif

      enddo


      te=ekine/xne          
      ti=ekini/xni

      call SERV_labelgen(mpe,node)
      call SERV_labelgen(n,label)
      
      
      if(n.eq.0) then 
        open(11,file=trim(data_out)//'/'//node//trim(fname),
     &     form='formatted')
c	 write(11) '# time  Te   Ti   (Ee-Ei)/(Ee+Ei)   Ee+Ei   Ne   Ni'
        write(11,1001) n*dt/wl
        write(11,1001) te
        write(11,1001) ti
        write(11,1001) (ekine-ekini)/(ekine+ekini)
        write(11,1001) ekine+ekini
        write(11,1001) xne
        write(11,1001) xni
      else 
      	open(11,file=trim(data_out)//'/'//node//trim(fname),
     & form='formatted', position='append')
        write(11,1001) n*dt/wl
        write(11,1001) te
        write(11,1001) ti
        write(11,1001) (ekine-ekini)/(ekine+ekini)
        write(11,1001) ekine+ekini
        write(11,1001) xne
        write(11,1001) xni
      endif
      close(11)


      if( n.eq.0.OR.n.eq.nprf) then    ! with pfield period 
        nprf=nprf+dnprf

        open(11,file=trim(data_out)
     &     //'/'//node//'e'//trim(fnames)//label,form='formatted')
c        write(11) '# Ek[MeV]   dN/dE   Ntot=', ntot
        do l=0,399 
	   write(11,1001) ekinl+l*dekin/399, enwe(l)
        enddo
        close(11)

        open(11,file=trim(data_out)
     &     //'/'//node//'i'//trim(fnames)//label,form='formatted')
c        write(11) '# Ek[MeV]   dN/dE   Ntot=', ntot
        do l=0,399 
	   write(11,1001) ikinl+l*dikin/399, enwi(l)
        enddo
        close(11)
      endif

      deallocate(enwe)
      deallocate(enwi)
 

 1000 FORMAT(1x,1(1x, 1i10)) 
 1001 FORMAT(1x,1(1x, 1e12.6)) 


      end subroutine OUT_energy
