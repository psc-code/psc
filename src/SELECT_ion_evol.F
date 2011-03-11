      subroutine SELECT_ion_evol

      use PIC_variables
      use VLA_variables
      use PSC_globals

      implicit none
      integer tm,nc

      character*5 label,node

      real(kind=4) :: xi,yi,zi
      real(kind=4) :: pxi,pyi,pzi
      real(kind=4) :: qni,mni,cni,lni,wni


      write(6,*) '*****************'
      write(6,*) 'IONS'
      write(6,*) '*****************'
      write(6,*) 'SELECT TIME:'
      read(5,*) tm


      nc=0
      call SERV_labelgen(tm,label)
      do pec=0,xnpe*ynpe*znpe-1
         call SERV_labelgen(pec,node)
         open(11,file=trim(data_out)//'/'//node//'ion'//label,
     &        access='sequential',form='unformatted')
         read(11,end=3001) cori,shift_z
         read(11,end=3001) i1n,i1x,i2n,i2x,i3n,i3x
         read(11,end=3001) np
         read(11,end=3001) nmax,xnpe,ynpe
         read(11,end=3001) nprparti,dnprparti,nistep
         read(11,end=3001) dt,dx,dy,dz
         read(11,end=3001) i0,e0,b0,j0,rho0,n0,cori,shift_z
         read(11,end=3001) wl,wp,ld
         read(11,end=3001) alpha,beta,eta
         read(11,end=3001) vos,vt,cc,qq,mm,tt,eps0
 3000    read(11,end=3001) xi,yi,zi,pxi,pyi,pzi,qni,mni,cni,lni,wni
         nc=nc+1
         goto 3000
 3001    close(11)
      enddo


      open(12,file=trim(data_processed)//'/'//'I'//label//'.data',
     &     access='sequential',form='formatted')

      write(12,1000) nc
      write(12,1000) i1n,i1x,i2n,i2x,i3n,i3x
      write(12,1000) tm,np
      write(12,1000) nmax,xnpe,ynpe
      write(12,1000) nprparti,dnprparti,nistep
      write(12,1001) dt,dx,dy,dz
      write(12,1001) i0,e0,b0,j0,rho0,n0,cori,shift_z
      write(12,1001) wl,wp,ld
      write(12,1001) alpha,beta,eta
      write(12,1001) vos,vt,cc,qq,mm,tt,eps0

      call SERV_labelgen(tm,label)
      do pec=0,xnpe*ynpe*znpe-1
         call SERV_labelgen(pec,node)
         open(11,file=trim(data_out)//'/'//node//'ion'//label,
     &        access='sequential',form='unformatted')
         read(11,end=2001) cori,shift_z
         read(11,end=2001) i1n,i1x,i2n,i2x,i3n,i3x
         read(11,end=2001) np
         read(11,end=2001) nmax,xnpe,ynpe
         read(11,end=2001) nprparti,dnprparti,nistep
         read(11,end=2001) dt,dx,dy,dz
         read(11,end=2001) i0,e0,b0,j0,rho0,n0,cori,shift_z
         read(11,end=2001) wl,wp,ld
         read(11,end=2001) alpha,beta,eta
         read(11,end=2001) vos,vt,cc,qq,mm,tt,eps0
 2000    read(11,end=2001) xi,yi,zi,pxi,pyi,pzi,qni,mni,cni,lni,wni
         write(12,1001) xi
         write(12,1001) yi
         write(12,1001) zi
         write(12,1001) pxi
         write(12,1001) pyi
         write(12,1001) pzi
         write(12,1001) qni
         write(12,1001) mni
         write(12,1001) cni
         write(12,1001) lni
         write(12,1001) wni
         goto 2000
 2001    close(11)
      enddo
      close(12)


 1000 FORMAT(1x, 1(1x, 1i10))
 1001 FORMAT(1x, 1(1x, 1e12.6))


      end subroutine SELECT_ion_evol
