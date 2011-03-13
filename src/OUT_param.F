
c DATA PROTOCOL ROUTINE

      subroutine OUT_param

      use VLA_variables
      use OUT_particles
      use OUT_fields
      use OUT_poyn
      use PSC_globals
      use PSC_domain
      use PSC_patch

      implicit none

      if (mpe.eq.0) then
         write(6,*) ' '
         write(6,*) 'BOX SIZE, PARTICLE NUMBER PER CELL'
         write(6,*) 'lengthx:',lengthx
         write(6,*) 'lengthy:',lengthy
         write(6,*) 'lengthz:',lengthz
         write(6,*) ' '
         write(6,*) 'BOUNDARY CONDITIONS'
         write(6,*) 'boundary_fld_x1:',boundary_fld_x1
         write(6,*) 'boundary_fld_x2:',boundary_fld_x2
         write(6,*) 'boundary_fld_y1:',boundary_fld_y1
         write(6,*) 'boundary_fld_y2:',boundary_fld_y2
         write(6,*) 'boundary_fld_z1:',boundary_fld_z1
         write(6,*) 'boundary_fld_z2:',boundary_fld_z2
         write(6,*) 'boundary_part_x:',boundary_part_x
         write(6,*) 'boundary_part_y:',boundary_part_y
         write(6,*) 'boundary_part_z:',boundary_part_z
         write(6,*) ' '
         write(6,*) 'PHYSICAL PARAMETERS'
         write(6,*) 'qq:',qq
         write(6,*) 'mm:',mm
         write(6,*) 'n0:',n0
         write(6,*) 'i0:',i0
         write(6,*) 'tt:',tt
         write(6,*) 'phi0:',phi0
         write(6,*) 'a0:',a0
         write(6,*) 'e0:',e0
         write(6,*) 'b0:',b0
         write(6,*) 'rho0:',rho0
         write(6,*) 'j0:',j0
         write(6,*) ' '
         write(6,*) 'NORMALIZATION PARAMETERS IN PHYSICAL UNITS'
         write(6,*) 'wp:',wp
         write(6,*) 'wl:',wl
         write(6,*) 'lw:',lw
         write(6,*) 'vt:',vt
         write(6,*) 'ld:',ld
         write(6,*) 'vos:',vos
         write(6,*) ' '
         write(6,*) 'NORMALIZATION PARAMETERS IN DIMENSIONLESS UNITS'
         write(6,*) 'alpha=wp/wl:',alpha
         write(6,*) 'beta=vt/c:',beta
         write(6,*) 'eta=vos/c:',eta
         write(6,*) ' '
         write(6,*) 'RESOLUTION IN PHYSICAL UNITS'
         write(6,*) 'dx:',dx*ld
         write(6,*) 'dy:',dy*ld
         write(6,*) 'dz:',dz*ld
         write(6,*) 'dt:',dt/wl
         write(6,*) ' '
         write(6,*) 'RESOLUTION IN DIMENSIONLESS UNITS'
         write(6,*) 'dx/ld:',dx
         write(6,*) 'dy/ld:',dy
         write(6,*) 'dz/ld:',dz
         write(6,*) 'wl*dt:',dt
         write(6,*) ' '
         write(6,*) 'JOB CONTROL PARAMETERS'
         write(6,*) 'nmax:',nmax
         write(6,*) 'nnp:',nnp
         write(6,*) 'np:',np
         write(6,*) 'tmnvf:',tmnvf
         write(6,*) 'tmxvf:',tmxvf
         write(6,*) 'tmnvp:',tmnvp
         write(6,*) 'tmxvp:',tmxvp
         write(6,*) 'nprf:',nprf
         write(6,*) 'dnprf:',dnprf
         write(6,*) 'nprparti:',nprparti
         write(6,*) 'dnprparti:',dnprparti
         write(6,*) 'nistep:',nistep
         write(6,*) 'xnpe: ',xnpe
         write(6,*) 'ynpe: ',ynpe
         write(6,*) 'znpe: ',znpe
         write(6,*) 'data_out:',trim(data_out)
         write(6,*) 'data_chk:',trim(data_chk)
         write(6,*) ' ' 
         write(6,*) 'Grid size output data'
         write(6,*) 'r1n:',r1n
         write(6,*) 'r1x:',r1x
         write(6,*) 'r2n:',r2n
         write(6,*) 'r2x:',r2x
         write(6,*) 'r3n:',r3n
         write(6,*) 'r3x:',r3x
         write(6,*) ' ' 
      endif

      end subroutine OUT_param
