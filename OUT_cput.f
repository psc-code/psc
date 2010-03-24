c WALL CLOCK TIME PROTOCOL PER TIME STEP

c cpu(1): WALL CLOCK time since program start
c cpu(2): WALL CLOCK time required for computation
c cpu(3): WALL CLOCK time required for communication
c cpu(4): WALL CLOCK time required for I/O
c cpu(5): LOCAL particle number


      if (mpe.eq.0) then

         write(6,*) ' '
         write(6,*) 'TIMESTEP:',n

         nitot=0
         do pec=0,npe-1
            cpu(1)=cpu_ary(1,pec)
            cpu(2)=cpu_ary(2,pec)
            cpu(3)=cpu_ary(3,pec)
            cpu(4)=cpu_ary(4,pec)
            cpu(5)=cpu_ary(5,pec)
            write(6,*) 'PE:',pec
            write(6,*) 'WT TOTAL:',sngl(cpu(1))
            write(6,*) 'WT COMPUTATION:',sngl(cpu(2))
            write(6,*) 'WT COMMUNICATION:',sngl(cpu(3))
            write(6,*) 'WT IN/OUT:',sngl(cpu(4))
            write(6,*) 'PNUM LOCAL:',sngl(cpu(5))
            nitot=nitot+cpu(5)
         enddo

         write(6,*) 'PNUM TOTAL:',nitot
         write(6,*) ' '

      endif
