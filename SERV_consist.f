c PARAMETER CONSISTENCY CHECK FOR THOSE PARAMETERS THAT WOULD
c CAUSE MASSIVE CODE FAILURE IF WRONGLY SETUP

      if (xnpe*ynpe*znpe.ne.npe) then
         if (mpe.eq.0) then
            write(6,*) 'NPE INCONSISTENT WITH DATA PARTITION!'
            write(6,*) 'REQUIRED NUMBER OF PEs:',xnpe*ynpe*znpe
            write(6,*) 'ALLOCATED NUMBER OF PEs:',npe
         endif
         call MPI_FINALIZE(info)
         stop 
      endif
      if (lengthx*lengthy*lengthz.eq.0.0) then
         if (mpe.eq.0) then
            write(6,*) 'ILLEGAL PHYSICAL BOX DIMENSIONS!'
         endif
         call MPI_FINALIZE(info)
         stop 
      endif
c      if (i1x-i1n+1.gt.i1tot) then
c         if (mpe.eq.0) then
c            write(6,*) 'ILLEGAL GRID. Enlarge i1tot!'
c         endif
c         call MPI_FINALIZE(info)
c         stop 
c      endif
c      if (i2x-i2n+1.gt.i2tot) then
c         if (mpe.eq.0) then
c            write(6,*) 'ILLEGAL GRID. Enlarge i2tot!'
c         endif
c         call MPI_FINALIZE(info)
c         stop 
c      endif
c      if (i3x-i3n+1.gt.i3tot) then
c         if (mpe.eq.0) then
c            write(6,*) 'ILLEGAL GRID. Enlarge i3tot!'
c         endif
c         call MPI_FINALIZE(info)
c         stop 
c      endif
