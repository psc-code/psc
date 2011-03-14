c CHECKPOINTING ROUTINE FOR SIMULATION DATA


      subroutine SERV_write(timestep, niloc, p_niloc,
     &     jxi, jyi, jzi, ex, ey, ez, hx, hy, hz)

      use VLA_variables, only: shift_c,shift_z
      use OUT_params
      use OUT_particles
      use OUT_fields
      use PIC_bin_collision ! for I/O params
      use OUT_poyn
      use PSC_globals
      use PSC_patch

      implicit none

      integer timestep
      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ex, ey, ez
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: hx, hy, hz

      integer i1,i2,i3
      integer i1a,i1e,i2a,i2e,i3a,i3e
      character*5 node

      call SERV_labelgen(mpe,node)

      open(10,file=trim(data_chk)//'/'//node//'CPNEW',
     &     access='sequential',form='unformatted')

      write(10) timestep
      write(10) shift_c,shift_z
      write(10) nprf,nprc
      write(10) nprparti
      write(10) tmnvf,tmxvf
      write(10) tmnvp,tmxvp

      call OUT_poyn_serv_write

      write(10) niloc,cori
      write(10) i1mn,i1mx,i2mn,i2mx,i3mn,i3mx


      do i1=0,11*niloc+10,100
         i1e=min(i1+99,11*niloc+10)
         write(10) (p_niloc(i1a),i1a=i1,i1e)
      enddo

      do i3=i3mn-rd3,i3mx+rd3,100
         i3e=min(i3+99,i3mx+rd3)
         do i2=i2mn-rd2,i2mx+rd2,100
            i2e=min(i2+99,i2mx+rd2)
             do i1=i1mn-rd1,i1mx+rd1,100
                i1e=min(i1+99,i1mx+rd1)
                write(10) (((jxi(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((jyi(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((jzi(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((ex(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((ey(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((ez(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((hx(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((hy(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((hz(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
            enddo
         enddo
      enddo

      call OUT_FIELDS_serv_write

      close(10)

      end subroutine SERV_write
