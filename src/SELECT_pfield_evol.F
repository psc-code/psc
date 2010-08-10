      subroutine SELECT_pfield_evol

      use PIC_variables
      use VLA_variables

      implicit none
      integer tm,index,ic
      integer :: i1e,i2e,i3e,i1a,i2a,i3a

      character*5 label,node
      character(6) fieldname
      character (len=6), dimension (34) :: f_name    ! changed by ab

      real(kind=4) :: x0,y0,z0,xr,yr,zr,rot,field
      real(kind=4),allocatable,dimension(:,:,:) :: t_serie


      f_name(1)='net'
      f_name(2)='nit'
      f_name(3)='nnt'
      f_name(4)='jxit'
      f_name(5)='jyit'
      f_name(6)='jzit'
      f_name(11)='ext'
      f_name(12)='eyt'
      f_name(13)='ezt'
      f_name(14)='bxt'
      f_name(15)='byt'
      f_name(16)='bzt'
      f_name(17)='hxt'                  ! ab
      f_name(18)='hyt'                  ! ab
      f_name(19)='hzt'                  ! ab
      f_name(20)='jxexit'
      f_name(21)='jyeyit'
      f_name(22)='jzezit'
      f_name(23)='poyxt'
      f_name(24)='poyyt'
      f_name(25)='poyzt'
      f_name(26)='ex2t'
      f_name(27)='ey2t'
      f_name(28)='ez2t'
      f_name(29)='bx2t'
      f_name(30)='by2t'
      f_name(31)='bz2t'
      f_name(32)='hx2t'                 ! ab
      f_name(33)='hy2t'                 ! ab
      f_name(34)='hz2t'                 ! ab


      write(6,*) 'SELECT FIELD NAME:'
      read(5,*) fieldname
      write(6,*) 'SELECT TIME:'
      read(5,*) tm


      allocate(t_serie(i1n:i1x,i2n:i2x,i3n:i3x))
      t_serie=0.0


      index=0
      do ic=1,34                         ! changed by ab
        if (fieldname.eq.f_name(ic)) then 
          index=ic
        endif
      enddo
      if (index.eq.0) then
        write(6,*) '********************'
        write(6,*) 'FIELD NOT AVAILABLE!'
        write(6,*) '********************'
        goto 999
      endif


      call SERV_labelgen(tm,label)
      do pec=0,xnpe*ynpe*znpe-1

         call SERV_labelgen(pec,node)
         open(11,file=trim(data_out)//'/'//node//'pfield'//label,
     &        access='sequential',form='unformatted')

         read(11) i1n,i1x,i2n,i2x,i3n,i3x
         read(11) np
         read(11) nmax,xnpe,ynpe
         read(11) nprf,dnprf
         read(11) dt,dx,dy,dz,shift_z
         read(11) i0,e0,b0,j0,rho0,n0
         read(11) wl,wp,ld
         read(11) alpha,beta,eta
         read(11) vos,vt,cc,qq,mm,tt,eps0
         read(11) i1mn,i1mx,i2mn,i2mx,i3mn,i3mx,shift_z
         read(11) r1n,r1x,r2n,r2x,r3n,r3x
         read(11) x0,y0,z0,rot

         do i3=i3mn,i3mx
            do i2=i2mn,i2mx
               do i1=i1mn,i1mx

                  xr=i1*dx
                  yr=cos(rot)*(i2*dy-y0)-sin(rot)*(i3*dz-z0)+y0
                  zr=cos(rot)*(i3*dz-z0)+sin(rot)*(i2*dy-y0)+z0

                  if (((r3n-1)*dz.le.zr.and.zr.le.(r3x+1)*dz).and.
     &                ((r2n-1)*dy.le.yr.and.yr.le.(r2x+1)*dy).and.
     &                ((r1n-1)*dx.le.xr.and.xr.le.(r1x+1)*dx)) then

                     do ic=1,index
                        read(11) field
                     enddo
                     t_serie(i1,i2,i3)=field
                     do ic=index+1,34            ! changed by ab
                        read(11) field
                     enddo

                  endif

               enddo
            enddo
         enddo

         close(11)

      enddo


      open(11,file=trim(data_processed)//'/'//trim(fieldname)//label
     & //'.data',
     &     access='sequential',form='formatted')
      open(12,file=trim(data_processed)//'/'//trim(fieldname)//label
     & //'.vtk',
     &     access='sequential',form='formatted')


      write(12,'(A26)') '# vtk DataFile Version 2.0'
      write(12,'(A24)') 'GridSetup generated file'
      write(12,'(A5)') 'ASCII'
      write(12,'(A25)') 'DATASET STRUCTURED_POINTS'
      write(12,'(A10,A1,I4,A1,I4,A1,I4)') 'DIMENSIONS',' ',
     &i1x-i1n+1,' ',i2x-i2n+1,' ',i3x-i3n+1
      write(12,'(A7,A1,E10.2,A1,E10.2,A1,E10.2)') 'SPACING',' ',
     &dx,' ',dy,' ',dz
      write(12,'(A6,A1,E10.2,A1,E10.2,A1,E10.2)') 'ORIGIN',' ',
     &(i1n*dx),' ',(i2n*dy),' ',(i3n*dz)
      write(12,'(A1)') ' '
      write(12,'(A10,A1,I8)') 'POINT_DATA',' ',
     &(i1x-i1n+1)*(i2x-i2n+1)*(i3x-i3n+1)
      write(12,'(A23)') 'SCALARS FIELDS double 1'
      write(12,'(A20)') 'LOOKUP_TABLE default'

      write(11,1000) i1n,i1x,i2n,i2x,i3n,i3x
      write(11,1000) tm,np
      write(11,1000) nmax,xnpe,ynpe
      write(11,1000) nprf,dnprf
      write(11,1001) dt,dx,dy,dz,shift_z
      write(11,1001) i0,e0,b0,j0,rho0,n0
      write(11,1001) wl,wp,ld
      write(11,1001) alpha,beta,eta
      write(11,1001) vos,vt,cc,qq,mm,tt,eps0

      do i3=i3n,i3x
         do i2=i2n,i2x
            do i1=i1n,i1x
               write(11,1001) t_serie(i1,i2,i3)
               write(12,1001) t_serie(i1,i2,i3)
            enddo
         enddo
      enddo

      close(11)
      close(12)

 999  continue


 1000 FORMAT(1x, 1(1x, 1i6))
 1001 FORMAT(1x, 1(1x, 1e12.5))


      end subroutine SELECT_pfield_evol
