
      module PSC_params

c i0: irradiance in Wum**2/m**2
c n0: peak atom density in m**3
c lw: laser wavelength
c wl: laser frequency
c qq: reference charge in As
c mm: reference mass in kg
c tt: atom temperature in J
c cc: light velocity in m/s
c eps0: eps0 in As/Vm
c wp: plasma frequency
c vt: thermal velocity
c ld: plasma wavelength

      real(kind=8) i0,n0,lw,wl
      real(kind=8) qq,mm,tt,cc,eps0
      real(kind=8) ld,wp,vos,vt
      real(kind=8) e0,b0,j0,rho0, phi0, a0

      real(kind=8) beta

      end module PSC_params
