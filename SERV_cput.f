c CPU TIME AND PARTICLE NUMBER BROADCAST


      cpu_ary(1,mpe)=cpuf-cpus
      cpu_ary(2,mpe)=cpucomp
      cpu_ary(3,mpe)=cpumess
      cpu_ary(4,mpe)=cpuinou
      cpu_ary(5,mpe)=niloc

      do pec=0,npe-1
         cpu(1)=cpu_ary(1,pec)
         cpu(2)=cpu_ary(2,pec)
         cpu(3)=cpu_ary(3,pec)
         cpu(4)=cpu_ary(4,pec)
         cpu(5)=cpu_ary(5,pec)
         call MPI_BCAST(cpu,5,MPI_DOUBLE_PRECISION,
     &                  pec,MPI_COMM_WORLD,info)
         cpu_ary(1,pec)=cpu(1)
         cpu_ary(2,pec)=cpu(2)
         cpu_ary(3,pec)=cpu(3)
         cpu_ary(4,pec)=cpu(4)
         cpu_ary(5,pec)=cpu(5)
      enddo


c LARGEST CPU-TIME FROM ALL NODES


      cpue=cpu_ary(1,0)
      do pec=0,npe-1
         if (cpu_ary(1,pec).gt.cpue) then
            cpue=cpu_ary(1,pec)
         endif
      enddo
