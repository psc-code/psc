      subroutine SELECT_count_evol

      use OUT_particles
      use VLA_variables
      use PSC_globals
      use PSC_params
      use PSC_domain
      use PSC_patch

      implicit none
      integer :: tm,nc
      integer :: count_in,count_ou,pec

      character*5 label,node

      real(kind=4) :: xi,yi,zi
      real(kind=4) :: pxi,pyi,pzi
      real(kind=4) :: qni,mni,cni,lni,wni
      real(kind=4) :: count_in_sum,count_ou_sum,time

      write(6,*) '*****************'
      write(6,*) 'COUNTING'
      write(6,*) '*****************'
      write(6,*) 'SELECT TIME:'
      read(5,*) tm


      nc=0
      count_in_sum=0
      count_ou_sum=0
      call SERV_labelgen(tm,label)
      do pec=0,xnpe*ynpe*znpe-1
         call SERV_labelgen(pec,node)
         open(11,file=trim(data_out)//'/'//node//'count'//label,
     &        access='sequential',form='unformatted')
 3000    read(11,end=3001) count_in,count_ou,time
         nc=nc+1
         count_in_sum=count_in_sum+count_in
         count_ou_sum=count_ou_sum+count_ou
         goto 3000
 3001    close(11)
      enddo


      open(12,file=trim(data_processed)//'/'//'COUNT'//label//'.data',
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
      write(12,1001) count_in_sum,count_ou_sum,time

      close(12)


 1000 FORMAT(1x, 1(1x, 1i10))
 1001 FORMAT(1x, 1(1x, 1e12.6))


      end subroutine SELECT_count_evol
