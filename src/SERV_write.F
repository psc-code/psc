c CHECKPOINTING ROUTINE FOR SIMULATION DATA


      subroutine SERV_write(timestep)

      use VLA_variables
      use PIC_variables
      use fields

      implicit none

      integer i1a,i1e,i2a,i2e,i3a,i3e
      integer timestep
      character*5 node


      call SERV_labelgen(mpe,node)

      open(10,file=trim(data_chk)//'/'//node//'CPNEW',
     &     access='sequential',form='unformatted')

      write(10) timestep+1
      write(10) shift_c,shift_z
      write(10) nprf,nprc
      write(10) nprparti
      write(10) tmnvf,tmxvf
      write(10) tmnvp,tmxvp
      write(10) tmnvc,tmxvc

      call OUT_poyn_serv_write

      write(10) niloc,nialloc,cori
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
                write(10) (((ne(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((ni(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
                write(10) (((nn(i1a,i2a,i3a),i1a=i1,i1e),
     &                       i2a=i2,i2e),i3a=i3,i3e)
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

#ifdef USE_PML
      do i3=i3mn,i3mx,100
         i3e=min(i3+99,i3mx)
         do i2=i2mn,i2mx,100
            i2e=min(i2+99,i2mx)
            do i1=i1mn,i1mx,100
               i1e=min(i1+99,i1mx)
               write(10) (((hxt(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((hyt(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((hzt(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((ex2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((ey2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((ez2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((bx2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((by2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((bz2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((hx2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((hy2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((hz2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((jxexit(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((jyeyit(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((jzezit(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((poyxt(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((poyyt(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
               write(10) (((poyzt(i1a,i2a,i3a),i1a=i1,i1e),
     &                      i2a=i2,i2e),i3a=i3,i3e)
            enddo
         enddo
      enddo
#endif

      close(10)

      end subroutine SERV_write
