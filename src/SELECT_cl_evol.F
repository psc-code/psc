#include "config.h"

      subroutine SELECT_cl_evol

      use VLA_variables
      use out_poyn
      use PSC_globals
      use PSC_patch

      implicit none
      integer tm,time,tstart,tend,tstep
      integer nodei,nodej,nodek

      character*5 label,node

      real(kind=8) total,absor,unit, glo_sum_abs,glo_abs      ! ab
      real(kind=8) glo_sum_enEXt,glo_sum_enEYt,glo_sum_enEZt
      real(kind=8) glo_sum_enBXt,glo_sum_enBYt,glo_sum_enBZt
      real(kind=8) glo_sum_enHXt,glo_sum_enHYt,glo_sum_enHZt  ! ab
      real(kind=8) glo_sum_ent,glo_sum_jet
      real(kind=8) glo_sum_poxt,glo_sum_poyt,glo_sum_pozt
      real(kind=8) glo_sum_poynit,glo_sum_poynot              ! ab

      real(kind=8) glo_fluxit
      real(kind=8) glo_fluxot
      real(kind=8) glo_enEXt
      real(kind=8) glo_enEYt
      real(kind=8) glo_enEZt
      real(kind=8) glo_enBXt
      real(kind=8) glo_enBYt
      real(kind=8) glo_enBZt
      real(kind=8) glo_enHXt          ! ab
      real(kind=8) glo_enHYt          ! ab
      real(kind=8) glo_enHZt          ! ab

      real(kind=8) glo_ent
      real(kind=8) glo_jet
      real(kind=8) glo_poxt
      real(kind=8) glo_poyt
      real(kind=8) glo_pozt

      real(kind=8) glo_poynit         ! ab
      real(kind=8) glo_poynot         ! ab

      real(kind=8) sum_enEXt
      real(kind=8) sum_enEYt
      real(kind=8) sum_enEZt
      real(kind=8) sum_enBXt
      real(kind=8) sum_enBYt
      real(kind=8) sum_enBZt
      real(kind=8) sum_enHXt          ! ab
      real(kind=8) sum_enHYt          ! ab
      real(kind=8) sum_enHZt          ! ab

      real(kind=8) sum_ent
      real(kind=8) sum_poxt
      real(kind=8) sum_poyt
      real(kind=8) sum_pozt
      real(kind=8) sum_jet

      real(kind=8) sum_poynit         ! ab
      real(kind=8) sum_poynot         ! ab


c PREPARATION


      write(6,*) ' '
      write(6,*) 'POYNTING LAW'
      write(6,*) ' '
      write(6,*) 'SELECT INITIAL TIME:'
      read(5,*) tstart
      write(6,*) 'SELECT FINAL TIME:'
      read(5,*) tend
      write(6,*) 'SELECT TIME STEP:'
      read(5,*) tstep


      open(11,file=trim(data_processed)//'/'//'poynting'//'.data',
     &     access='sequential',form='formatted')

      write(11,*) ' '
      write(11,*) 'NODE, GRID SIZE, BOX SIZE'
      write(11,*) 'i1n:',i1n
      write(11,*) 'i1x:',i1x
      write(11,*) 'i2n:',i2n
      write(11,*) 'i2x:',i2x
      write(11,*) 'i3n:',i3n
      write(11,*) 'i3x:',i3x
      write(11,*) ' '
      write(11,*) 'JOB CONTROL PARAMETERS'
      write(11,*) 'nmax:',nmax
      write(11,*) 'np:',np
      write(11,*) 'nprf:',nprf
      write(11,*) 'dnprf:',dnprf
      write(11,*) 'xnpe: ',xnpe
      write(11,*) 'ynpe: ',ynpe
      write(11,*) ' '
      write(11,*) 'PHYSICAL PARAMETERS'
      write(11,*) 'qq:',qq
      write(11,*) 'mm:',mm
      write(11,*) 'n0:',n0
      write(11,*) 'i0:',i0
      write(11,*) 'tt:',tt
      write(11,*) 'e0:',e0
      write(11,*) 'b0:',b0
      write(11,*) 'rho0:',rho0
      write(11,*) 'j0:',j0
      write(11,*) ' '
      write(11,*) 'NORMALIZATION PARAMETERS IN PHYSICAL UNITS'
      write(11,*) 'wp:',wp
      write(11,*) 'wl:',wl
      write(11,*) 'vt:',vt
      write(11,*) 'ld:',ld
      write(11,*) 'vos:',vos
      write(11,*) ' '
      write(11,*) 'NORMALIZATION PARAMETERS IN DIMENSIONLESS UNITS'
      write(11,*) 'alpha=wp/wl:',alpha
      write(11,*) 'beta=vt/c:',beta
      write(11,*) 'eta=vos/c:',eta
      write(11,*) ' '
      write(11,*) 'RESOLUTION IN PHYSICAL UNITS'
      write(11,*) 'dx:',dx*ld
      write(11,*) 'dy:',dy*ld
      write(11,*) 'dz:',dz*ld
      write(11,*) 'dt:',dt/wl
      write(11,*) ' '
      write(11,*) 'RESOLUTION IN DIMENSIONLESS UNITS'
      write(11,*) 'dx/ld:',dx
      write(11,*) 'dy/ld:',dy
      write(11,*) 'dz/ld:',dz
      write(11,*) 'wl*dt:',dt
      write(11,*) ' ' 


      unit=1.0e3*eps0*e0*e0*ld*ld*ld
      glo_sum_enEXt=0.0   
      glo_sum_enEYt=0.0
      glo_sum_enEZt=0.0
      glo_sum_enBXt=0.0
      glo_sum_enBYt=0.0
      glo_sum_enBZt=0.0
      glo_sum_enHXt=0.0      ! ab
      glo_sum_enHYt=0.0      ! ab
      glo_sum_enHZt=0.0      ! ab
      glo_sum_ent=0.0
      glo_sum_poxt=0.0
      glo_sum_poyt=0.0
      glo_sum_pozt=0.0
      glo_sum_jet=0.0

      glo_sum_poynit=0.0    ! ab
      glo_sum_poynot=0.0    ! ab

      do tm=tstart,tend,tstep
         call SERV_labelgen(tm,label)

         glo_fluxit=0.0
         glo_fluxot=0.0
         glo_enEXt=0.0
         glo_enEYt=0.0
         glo_enEZt=0.0
         glo_enBXt=0.0
         glo_enBYt=0.0
         glo_enBZt=0.0
         glo_enHXt=0.0   ! ab
         glo_enHYt=0.0   ! ab
         glo_enHZt=0.0   ! ab
         glo_ent=0.0
         glo_poxt=0.0
         glo_poyt=0.0
         glo_pozt=0.0
         glo_jet=0.0

         glo_poynit=0.0
         glo_poynot=0.0

         do pec=0,xnpe*ynpe*znpe-1

            nodei=seg_i1(pec)
            nodej=seg_i2(pec) 
            nodek=seg_i3(pec) 

            call SERV_labelgen(pec,node)

            open(12,file=trim(data_out)//'/'//node//'poynting'//label,
     &           access='sequential',form='unformatted')
            read(12) time
            read(12) i1mn
            read(12) i1mx
            read(12) i2mn
            read(12) i2mx
            read(12) i3mn
            read(12) i3mx
            read(12) fluxit
            read(12) fluxot
            read(12) enEXt
            read(12) enEYt
            read(12) enEZt
            read(12) enHXt
            read(12) enHYt
            read(12) enHZt
            read(12) enBXt       ! ab
            read(12) enBYt       ! ab
            read(12) enBZt       ! ab
            read(12) ent
            read(12) poxt
            read(12) poyt
            read(12) pozt
            read(12) jet
            read(12) poynit      ! ab
            read(12) poynot      ! ab
            close(12)

            if (nodek==1) then
               glo_fluxit=glo_fluxit+fluxit
               glo_fluxot=glo_fluxot+fluxot
            endif

            glo_enEXt=glo_enEXt+enEXt
            glo_enEYt=glo_enEYt+enEYt
            glo_enEZt=glo_enEZt+enEZt
            glo_enHXt=glo_enHXt+enHXt
            glo_enHYt=glo_enHYt+enHYt
            glo_enHZt=glo_enHZt+enHZt
            glo_enBXt=glo_enBXt+enBXt
            glo_enBYt=glo_enBYt+enBYt
            glo_enBZt=glo_enBZt+enBZt

            glo_ent=glo_ent+ent
            glo_poxt=glo_poxt+poxt
            glo_poyt=glo_poyt+poyt
            glo_pozt=glo_pozt+pozt
            glo_jet=glo_jet+jet

            glo_poynit=glo_poynit+poynit    ! ab
            glo_poynot=glo_poynot+poynot    ! ab
         enddo

         glo_sum_enEXt=glo_sum_enEXt+glo_enEXt
         glo_sum_enEYt=glo_sum_enEYt+glo_enEYt
         glo_sum_enEZt=glo_sum_enEZt+glo_enEZt
         glo_sum_enBXt=glo_sum_enBXt+glo_enBXt
         glo_sum_enBYt=glo_sum_enBYt+glo_enBYt
         glo_sum_enBZt=glo_sum_enBZt+glo_enBZt
         glo_sum_enHXt=glo_sum_enHXt+glo_enHXt     ! ab
         glo_sum_enHYt=glo_sum_enHYt+glo_enHYt     ! ab
         glo_sum_enHZt=glo_sum_enHZt+glo_enHZt     ! ab

         glo_sum_ent=glo_sum_ent+glo_ent
         glo_sum_poxt=glo_sum_poxt+glo_poxt
         glo_sum_poyt=glo_sum_poyt+glo_poyt
         glo_sum_pozt=glo_sum_pozt+glo_pozt
         glo_sum_jet=glo_sum_jet+glo_jet

         glo_sum_poynit=glo_sum_poynit+glo_poynit   ! ab
         glo_sum_poynot=glo_sum_poynot+glo_poynot   ! ab

         total=glo_ent+glo_poxt+glo_poyt+glo_pozt+glo_jet  ! - changed

         if (glo_fluxit.gt.1.0e-20) then
            absor=1.0e2*(glo_fluxit-glo_fluxot)/glo_fluxit
         else
            absor=1.0e30
         endif

         glo_sum_abs = 1.0e2*(1-glo_sum_poynot/glo_sum_poynit)
         glo_abs = 1.0e2*(1-glo_poynot/glo_poynit)

         write(11,*) ' '
         write(11,*) 'GLOBAL POYNTING LAW'
         write(11,*) 'TIMESTEP=',time
         write(11,*) 'TIME=',1.0e15*time*dt/wl,'fs'
         write(11,*) '------------------------------------------'
         write(11,*) 'FLUXI at z=',i3n*dz*ld,':',glo_fluxit*unit,'mJ'
         write(11,*) 'FLUXO at z=',i3n*dz*ld,':',glo_fluxot*unit,'mJ'
         write(11,*) 'ABSOR IN % OF IRRAD:',absor
         write(11,*) '------------------------------------------'
         write(11,*) 'ENERGY DISTRIBUTION (LIGHT)'
         write(11,*) 'ENEX:',glo_enEXt*unit,'mJ'
         write(11,*) 'ENEY:',glo_enEYt*unit,'mJ'
         write(11,*) 'ENEZ:',glo_enEZt*unit,'mJ'
         write(11,*) 'ENBX:',glo_enBXt*unit,'mJ'
         write(11,*) 'ENBY:',glo_enBYt*unit,'mJ'
         write(11,*) 'ENBZ:',glo_enBZt*unit,'mJ'
         write(11,*) 'ENHX:',glo_enHXt*unit,'mJ'     ! ab
         write(11,*) 'ENHY:',glo_enHYt*unit,'mJ'     ! ab
         write(11,*) 'ENHZ:',glo_enHZt*unit,'mJ'     ! ab
         write(11,*) '------------------------------------------'
         write(11,*) 'TOTAL ENERGY DISTRIBUTION (LIGHT)'
         write(11,*) 'TOT ENEX:',glo_sum_enEXt*unit,'mJ'
         write(11,*) 'TOT ENEY:',glo_sum_enEYt*unit,'mJ'
         write(11,*) 'TOT ENEZ:',glo_sum_enEZt*unit,'mJ'
         write(11,*) 'TOT ENBX:',glo_sum_enBXt*unit,'mJ'
         write(11,*) 'TOT ENBY:',glo_sum_enBYt*unit,'mJ'
         write(11,*) 'TOT ENBZ:',glo_sum_enBZt*unit,'mJ'
         write(11,*) 'TOT ENHX:',glo_sum_enHXt*unit,'mJ'  ! ab
         write(11,*) 'TOT ENHY:',glo_sum_enHYt*unit,'mJ'  ! ab
         write(11,*) 'TOT ENHZ:',glo_sum_enHZt*unit,'mJ'  ! ab
         write(11,*) '------------------------------------------'
         write(11,*) 'ENERGY CONSERVATION (LIGHT)'
         write(11,*) '1. ENDEN:',glo_ent*unit,'mJ'
         write(11,*) '2. ECURX:',glo_poxt*unit,'mJ'
         write(11,*) '3. ECURY:',glo_poyt*unit,'mJ'
         write(11,*) '4. ECURZ:',glo_pozt*unit,'mJ'
         write(11,*) '5. JE:',glo_jet*unit,'mJ'
         write(11,*) '1 + 2 + 3 + 4 + 5 :',total*unit,'mJ'
         write(11,*) '------------------------------------------'
         write(11,*) 'TOT ENDEN:',glo_sum_ent*unit,'mJ'
         write(11,*) 'TOT ECURX:',glo_sum_poxt*unit,'mJ'
         write(11,*) 'TOT ECURY:',glo_sum_poyt*unit,'mJ'
         write(11,*) 'TOT ECURZ:',glo_sum_pozt*unit,'mJ'
         write(11,*) 'TOT JE:',glo_sum_jet*unit,'mJ'
         write(11,*) '------------------------------------------'
         write(11,*) 'GLOB SUM ABSORPTION IN %'                     ! ab
         write(11,*) 'POYN IN', glo_sum_poynit                      ! ab 
         write(11,*) 'POYN OUT', glo_sum_poynot                     ! ab 
         write(11,*) 'ABS' , glo_sum_abs                            ! ab
         write(11,*) '------------------------------------------'
         write(11,*) 'GLOB ABSORPTION IN %'                     ! ab
         write(11,*) 'POYN IN', glo_poynit                      ! ab 
         write(11,*) 'POYN OUT', glo_poynot                     ! ab 
         write(11,*) 'ABS' , glo_abs                            ! ab
         write(11,*) '------------------------------------------'
         write(11,*) ' '
      enddo


      do pec=0,xnpe*ynpe*znpe-1
         call SERV_labelgen(pec,node)

         sum_enEXt=0.0
         sum_enEYt=0.0
         sum_enEZt=0.0
         sum_enBXt=0.0
         sum_enBYt=0.0
         sum_enBZt=0.0
         sum_enHXt=0.0  ! ab
         sum_enHYt=0.0  ! ab
         sum_enHZt=0.0  ! ab
         sum_ent=0.0
         sum_poxt=0.0
         sum_poyt=0.0
         sum_pozt=0.0
         sum_jet=0.0

         sum_poynit=0.0   ! ab
         sum_poynot=0.0   ! ab

         nodei=seg_i1(pec)
         nodej=seg_i2(pec) 
         nodek=seg_i3(pec) 

         do tm=tstart,tend,tstep
            call SERV_labelgen(tm,label)

            open(12,file=trim(data_out)//'/'//node//'poynting'//label,
     &           access='sequential',form='unformatted')
            read(12) time
            read(12) i1mn
            read(12) i1mx
            read(12) i2mn
            read(12) i2mx
            read(12) i3mn
            read(12) i3mx
            read(12) fluxit
            read(12) fluxot
            read(12) enEXt
            read(12) enEYt
            read(12) enEZt
            read(12) enHXt
            read(12) enHYt
            read(12) enHZt
            read(12) enBXt   ! ab
            read(12) enBYt   ! ab
            read(12) enBZt   ! ab
            read(12) ent
            read(12) poxt
            read(12) poyt
            read(12) pozt
            read(12) jet
            read(12) poynit      ! ab
            read(12) poynot      ! ab
            close(12)

            total=ent+poxt+poyt+pozt+jet   ! - changed

            sum_enEXt=sum_enEXt+enEXt
            sum_enEYt=sum_enEYt+enEYt
            sum_enEZt=sum_enEZt+enEZt
            sum_enHXt=sum_enHXt+enHXt
            sum_enHYt=sum_enHYt+enHYt
            sum_enHZt=sum_enHZt+enHZt
            sum_enBXt=sum_enBXt+enBXt    ! ab
            sum_enBYt=sum_enBYt+enBYt    ! ab 
            sum_enBZt=sum_enBZt+enBZt    ! ab

            sum_ent=sum_ent+ent
            sum_poxt=sum_poxt+poxt
            sum_poyt=sum_poyt+poyt
            sum_pozt=sum_pozt+pozt
            sum_jet=sum_jet+jet

            sum_poynit = sum_poynit+poynit     ! ab
            sum_poynot = sum_poynot+poynot     ! ab

            write(11,*) ' '
            write(11,*) 'LOCAL POYNTING LAW'
            write(11,*) 'TIMESTEP=',time
            write(11,*) 'TIME=',1.0e15*time*dt/wl,'fs'
            write(11,*) 'NODEI=',nodei
            write(11,*) 'NODEJ=',nodej
            write(11,*) 'NODEK=',nodek
            write(11,*) '------------------------------------------'
            write(11,*) 'FLUXI at z=',i3mn*dz*ld,':',fluxit*unit,'mJ'
            write(11,*) 'FLUXO at z=',i3mn*dz*ld,':',fluxot*unit,'mJ'
            write(11,*) '------------------------------------------'
            write(11,*) 'ENERGY DISTRIBUTION (LIGHT)'
            write(11,*) 'ENEX:',enEXt*unit,'mJ'
            write(11,*) 'ENEY:',enEYt*unit,'mJ'
            write(11,*) 'ENEZ:',enEZt*unit,'mJ'
            write(11,*) 'ENHX:',enHXt*unit,'mJ'
            write(11,*) 'ENHY:',enHYt*unit,'mJ'
            write(11,*) 'ENHZ:',enHZt*unit,'mJ'
            write(11,*) 'ENBX:',enBXt*unit,'mJ'   ! ab
            write(11,*) 'ENBY:',enBYt*unit,'mJ'   ! ab
            write(11,*) 'ENBZ:',enBZt*unit,'mJ'   ! ab
            write(11,*) '------------------------------------------'
            write(11,*) 'TOTAL ENERGY DISTRIBUTION (LIGHT)'
            write(11,*) 'TOT ENEX:',sum_enEXt*unit,'mJ'
            write(11,*) 'TOT ENEY:',sum_enEYt*unit,'mJ'
            write(11,*) 'TOT ENEZ:',sum_enEZt*unit,'mJ'
            write(11,*) 'TOT ENBX:',sum_enBXt*unit,'mJ'
            write(11,*) 'TOT ENBY:',sum_enBYt*unit,'mJ'
            write(11,*) 'TOT ENBZ:',sum_enBZt*unit,'mJ'
            write(11,*) 'TOT ENHX:',sum_enHXt*unit,'mJ'  ! ab
            write(11,*) 'TOT ENHY:',sum_enHYt*unit,'mJ'  ! ab
            write(11,*) 'TOT ENHZ:',sum_enHZt*unit,'mJ'  ! ab
            write(11,*) '------------------------------------------'
            write(11,*) 'ENERGY CONSERVATION (LIGHT)'
            write(11,*) '1. ENDEN:',ent*unit,'mJ'
            write(11,*) '2. ECURX:',poxt*unit,'mJ'
            write(11,*) '3. ECURY:',poyt*unit,'mJ'
            write(11,*) '4. ECURZ:',pozt*unit,'mJ'
            write(11,*) '5. JE:',jet*unit,'mJ'
            write(11,*) '1 + 2 + 3 + 4 + 5 :',total*unit,'mJ'
            write(11,*) '------------------------------------------'
            write(11,*) 'TOT ENDEN:',sum_ent*unit,'mJ'
            write(11,*) 'TOT ECURX:',sum_poxt*unit,'mJ'
            write(11,*) 'TOT ECURY:',sum_poyt*unit,'mJ'
            write(11,*) 'TOT ECURZ:',sum_pozt*unit,'mJ'
            write(11,*) 'TOT JE:',sum_jet*unit,'mJ'
            write(11,*) '------------------------------------------'  
            write(11,*) 'LOCAL POYN IN', sum_poynit                  ! ab
            write(11,*) 'LOCAL POYN OUT', sum_poynot                 ! ab
            write(11,*) '------------------------------------------' ! ab
            write(11,*) ' '

         enddo
      enddo
      close(11)

      deallocate(seg_i1,seg_i2)

      return
      end
