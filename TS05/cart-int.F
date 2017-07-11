      subroutine cart1d(allvar,points,nvar,nvarp3,npts,
     $     allcart,intfirst)
	include  'param.inc'
	include 'intco1d.inc'
	include 'meter.inc'
        logical intfirst
	real points(3,npts)
c     in the following changed nk -> nkp1 - bkress
*       real allvar(nvarp3,ni,njp2,nk)
        real allvar(nvarp3,ni,njp2,nkp1)
        real allcart(nvar,npts),arot(3,3)
	data re/6.38e8/
C     
c
      print *,'intfirst=',intfirst

        if ( intfirst ) then
	pi = atan(1.)*4.
c
	do 827 j=1,njp1
	do 827 i=1,nip1
	   x2(i,j) = x2(i,j)/re
	   y2(i,j) = y2(i,j)/re
  827	continue
*	write (6,*) (x2(i,1),i=1,nip1)
*	write (6,*) (phi(i),i=1,nkp1)
	  do 808 j=1,nj
	  do 808 i=1,ni
	   x2mid(i,j+1) = 0.25*(x2(i,j)+x2(i+1,j)+x2(i,j+1)+
     $		x2(i+1,j+1))
	   y2mid(i,j+1) = 0.25*(y2(i,j)+y2(i+1,j)+y2(i,j+1)+
     $		y2(i+1,j+1))
  808	  continue
	 do 807 i=1,ni
	    x2mid(i,1) = x2mid(i,2)
	    x2mid(i,njp2) = x2mid(i,njp1)
	    y2mid(i,1) = 0.
	    y2mid(i,njp2) = 0.
  807	continue
        do j=1,njp2
           psimid_in(j) = atan2(y2mid(1,j),x2mid(1,j))
           psimid_out(j) = atan2(y2mid(ni,j),x2mid(ni,j))
        enddo
	  do 809 k=1,nk
	    phimid(k) = 0.5*(phi(k)+phi(k+1))
  809	  continue
	    phimid(nkp1) = phimid(1) + 2.*pi
	  do 810 k=1,nk
	  do 810 j=1,nj
	  do 810 i=1,ni
	    x3mid(i,j,k+1) = x2mid(i,j)
	    y3mid(i,j,k+1) = y2mid(i,j)*cos(phimid(k))
	    z3mid(i,j,k+1) = y2mid(i,j)*sin(phimid(k))
  810	  continue
	do 812 j=1,njp2
	do 812 i=1,ni
	   x3mid(i,j,1) = x3mid(i,j,nkp1)
	   x3mid(i,j,nkp2) = x3mid(i,j,2)
	   y3mid(i,j,1) = y3mid(i,j,nkp1)
	   y3mid(i,j,nkp2) = y3mid(i,j,2)
	   z3mid(i,j,1) = z3mid(i,j,nkp1)
	   z3mid(i,j,nkp2) = z3mid(i,j,2)
  812	continue
        npt2 = npts

          print *,'calling intco1d'

	  call intco1d(points,npt2)
*	write (6,*) 'O.K. after coefficients, Biil'
c
        intfirst = .false.
c
        endif
C
C       PROCESSING THE INDIVIDUAL DATA
C
	nt = ni*njp2*nkp1
	write (6,*) 'interpolating data at this time step'
        call interp1d(allcart,allvar,ni,njp2,nkp1,nt,npts,nvar)
C
C
C
C
        return
	END






