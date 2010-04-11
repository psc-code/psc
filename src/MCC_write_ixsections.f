c     WRITE IONIZATION X-SECTIONS INTO A FILE FOR CONTROL PURPOSES ONLY

      subroutine MCC_write_ixsections

      use MCC_variables
      
      implicit none

      real(kind=8) :: e0,e1,f,ek
      real(kind=8) :: sigma,sv
      integer      :: i,j,k,mat
      integer      :: NPOINTS

      NPOINTS=100

      do mat=1,NMAT

         write(*,*) ' Electron Impact Ionization X-Section for : ',
     &matname(mat)
         write(*,*) ' '            
         do i=0,n_xstable(mat)-1
            write(*,*) ' Charge state ',i
            write(*,*) ' Ekin[eV]  sigma[m^2]  sigma v[m^3/s]'
            e0=xstable_t(mat,i)
            e1=xstable_e(mat,i,xstable_n(mat,i)-1)
            f=(e1/e0)**(1.0/NPOINTS)
            
            ek=e0
            do j=0,NPOINTS
               call MCC_ixsection(mat,i,ek,sigma,sv)
               write(*,*) ek,sigma, sv
               ek=ek*f
            enddo
            write(*,*) ' '
         enddo
      enddo

      end subroutine MCC_write_ixsections





