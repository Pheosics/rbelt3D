c boundx.inc
      PARAMETER (ntable=20000)
      REAL njp,vjx,vjy,vjz,cs,bjx,bjy,bjz,bjt,vsw,deltat
      INTEGER index
      LOGICAL table, sm
      COMMON /WIND/
     $     mjin,mjinp1,mjinp2,mjwall,mjwp1,mjwp2,
     $     tzero,
     $     rho_0,vx_0,vy_0,vz_0,c_0,bx_0,by_0,bz_0,
     $     frmat(3,3),
     $	   buff_in(2,nj,nk,9),
     $     buff_b(njp1,nkp1,9),buff_bx(njp1,nk),
     $	   polegse(3),tilt_rot(3,3),tilt_angle,
     $	   table,itable,index,deltat,
     $	   njp(ntable), vjx(ntable),vjy(ntable),vjz(ntable),
     $     ut(ntable),cs(ntable),bjX(ntable),bjy(ntable),
     $     bjz(ntable),bjt(ntable),zangle(ntable),
     $     by_coef,bz_coef,bx_zero,sm
	common/buff2d/ buff_h(2,nj,9,2)
!$OMP THREADPRIVATE (/buff2d/)
c
c	mjin = the number of inflowing j-faces.
c       mjinp1 = the j-value of the front corner?
c	mjwall = the j-value of the rear corner
