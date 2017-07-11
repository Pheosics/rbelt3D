c     const.inc - rbelt constants

      real Re,m0,c,charge,q,ntg,vmsvcm,wrest,ec2sc,charge_sign
      real tfactor,ffactor,b0,raddeg,pi,svv,rgdt2gv

c     scaling constants (http://physics.nist.gov/)
c     define charge state & mass here ********************************
      parameter(charge = 1.)		! elementary charge (an integer)
      parameter(charge_sign = 1.)    ! -1. for time reversal or electrons
      parameter(m0=1.67262171e-24)		!proton rest mass, g
*      parameter(m0=9.1093826e-28)  !electron rest mass, g
c     ****************************************************************
      parameter(Re=6378.137e5)		! earth eq. radius, in cm (IUGG)
c     note, polar radius is ~ 6356.75e5 cm
c     mean radius (r=dsqrt(A/4*pi)) is ~ 6371.3e5 cm
c     quadratic mean radius (q=dsqrt((3r_e^2+r_p^2)/4) is ~ 6372.795e5 cm
      parameter(c=2.99792458e10)		! speed of light, cm/s
      parameter(ec2sc=3.e9*1.60217653e-19) ! elem. charge to stat coul.
      parameter(q=charge*ec2sc)		! charge, in esu
      parameter(ffactor=charge_sign*q*Re/(m0*c*c))
      parameter(tfactor=c/Re)

c     additional useful constants
      parameter(wrest=m0*c*c*3.e9/1.e13/ec2sc) !rest energy, MeV
      parameter(ntg=1.e-5)			!converts nT to gauss
      parameter(vmsvcm=1.e6/c)			!converts V/m to sV/cm
      parameter(svv=c/1.e8)			!sV to V
c     dimensionless ridigity to GV
      parameter(rgdt2gv=svv*Re/(q*Re/(m0*c*c))/1.e9)
      parameter(pi = 3.14159265358979)
      parameter(raddeg = 180./pi)

c     subtracted dipole field (added to pert. fields in linterp)
c     ****************************************************************
c     must be the same as b0 used in lfm_sample.
c     For best results, should also be the same as
c     GEOQMU in INPUT1 file used for lfm run
*      parameter(b0 = .27500*ffactor)
      parameter(b0 = .30000*ffactor)
*      parameter(b0 = .31100*ffactor)
c     ****************************************************************
