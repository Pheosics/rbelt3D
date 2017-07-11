
	subroutine intco1d(points,npts)
	
	include 'param.inc'
      dimension arot(3,3)
	parameter (ninj=ni*nj, ninjnk=ni*nj*nk, ninjp2=ni*njp2)
	include 'meter.inc'
	include 'intco1d.inc'
	
	
      real points(3,npts)
*	common /scrch1/ dxiq(2,nim1,njp2),dxjq(2,ni,njp1)
c	3rd index in dxjq exceeds njp1 below! changed to njp2 - bkress
	common /scrch1/ dxiq(2,nim1,njp2),dxjq(2,ni,njp2)
	real xdum(2)
	real eta1(num1d),chi1(num1d)
	real xi1(num1d),pv(3)
      real newdif
c	arrays below not used except phicart - bkress
*	real phicart(num1d),ycart1(num1d),distsq(num1d),phidif(nk)
	real phicart(num1d)
c	additional var. are used - bkress
*	real pi,olddif,xzero,yzero,psicart,rzero,psx
*	real dp,xouter,youter,router,si,sj,ssi,ssj
*	real dx,dy,xchi,xetachi,yeta,ychi,yetachi
*	real aa,bb,cc,disc,bt,etainc,chiinc
*	integer iint,ntries,jint,iold,jold,ii,jj,k,jinner,jouter

*	common /all/arot,xdum,eta1,chi1,xi1,pv,newdif,phicart
*	common /all2/pi,olddif,xzero,yzero,psicart,rzero,psx,jouter
*	common /all3/dp,xouter,youter,router,jinner,si,sj,ssi,ssj
*	common /all4/dx,dy,xchi,xetachi,yeta,ychi,yetachi
*	common /all5/aa,bb,cc,disc,bt,etainc,chiinc
*	common /all6/iint,ntries,jint,iold,jold,ii,jj,k
c
c       bkress (temp. hack)
        points(1,515140)=-2.3
c
	pi = atan(1.0)*4.0
c
	do 10 k=1,npts
   	phicart(k) = atan2(points(3,k),points(2,k))
   10	if (phicart(k).lt. phimid(1)) phicart(k) = 2.*pi+phicart(k)
c
c
	do 30 k=1,npts
c
           olddif = phicart(k) - phimid(1)
	do 20 n=1,nk
           newdif = phicart(k) - phimid(n+1)
           if ( newdif .le. 0.0 .and. olddif .gt. 0.0 ) then
c            here's the cell
             xi1(k) = float( mod(n+nk-1,nk) + 1) + 
     $             olddif/(olddif-newdif)
              go to 21
           endif
           olddif = newdif
 20     continue
 21     continue
c
c
   30	continue
c
c
	do 200 k=1,njp2
	do 200 j=1,nim1
	  dxiq(1,j,k) = x2mid(j+1,k)-x2mid(j,k)
	  dxiq(2,j,k) = y2mid(j+1,k)-y2mid(j,k)
  200	continue
	do 201 k=1,njp1
	do 201 j=1,ni
	  dxjq(1,j,k) = x2mid(j,k+1)-x2mid(j,k)
	  dxjq(2,j,k) = y2mid(j,k+1)-y2mid(j,k)
  201	continue

c      temporary bad fix until we figure out what was intended - bkress
       dxjq(1,1,26)=dxjq(1,1,25)
       dxjq(2,1,26)=dxjq(2,1,25)

c
c
	iint = 1
	do 400 k=1,npts
c
	xzero = points(1,k)
	yzero = sqrt(points(2,k)**2+points(3,k)**2)
c
c
        psicart = atan2(yzero,xzero)
        rzero = xzero**2 + yzero**2
c
c   test for outer bound
c
           olddif = psicart - psimid_out(1)
	do 70 n=1,njp1
           newdif = psicart - psimid_out(n+1)
           if ( newdif .le. 0.0 .and. olddif .ge. 0.0 ) then
c            here's the cell
              psx = float(n) + 
     $             olddif/(olddif-newdif)
              go to 71
           endif
           olddif = newdif
 70     continue
 71     continue
        jouter = psx
        dp = psx - float(jouter)
        xouter = x2mid(ni,jouter)*(1.-dp) + x2mid(ni,jouter+1)*dp
        youter = y2mid(ni,jouter)*(1.-dp) + y2mid(ni,jouter+1)*dp
        router = xouter**2 + youter**2
cccc
c   test for inner bound
c
           olddif = psicart - psimid_in(1)
	do 72 n=1,njp1
           newdif = psicart - psimid_in(n+1)
           if ( newdif .le. 0.0 .and. olddif .ge. 0.0 ) then
c            here's the cell
              psx = float(n) + 
     $             olddif/(olddif-newdif)
              go to 73
           endif
           olddif = newdif
 72     continue
 73     continue
        jinner = psx
        dp = psx - float(jinner)
        xinner = x2mid(1,jinner)*(1.-dp) + x2mid(1,jinner+1)*dp
        yinner = y2mid(1,jinner)*(1.-dp) + y2mid(1,jinner+1)*dp
        rinner = xinner**2 + yinner**2
cccc
ccc
cc
c       check to see if its in bounds
cc
ccc
cccc
	IF ( rzero .gt. rinner .and. rzero .lt. router) then
c
c
	ntries = 1
	iint=1
	jint = jinner
  417	continue
  
c
	do 399 kounter=1,200
	  iold = iint
	  jold = jint
	  xdum(1) = xzero - x2mid(iint,jint)
	  xdum(2) = yzero - y2mid(iint,jint)
	  si = sign(0.5,cross(dxiq(1,iint,jint),xdum))
	  sj = sign(0.5,cross(xdum,dxjq(1,iint,jint)))
	  if ( si + sj .gt. 0.75 ) then
	     xdum(1) = xzero - x2mid(iint+1,jint+1)
	     xdum(2) = yzero - y2mid(iint+1,jint+1)
	     ssi = sign(0.5,cross(xdum,dxiq(1,iint,jint+1)))
	     ssj = sign(0.5,cross(dxjq(1,iint+1,jint),xdum))
	     if ( ssi + ssj .gt. 0.75 ) then
		eta1(k) = float(iint)+0.5
		chi1(k) = float(jint)+0.5
	 	go to 1177   ! found cell !! ##################################   
	     else
		if ( ssi .lt. 0. ) jint = max0(1,min0(njp1,jint+1))
		if ( ssj .lt. 0. ) iint = max0(1,min0(nim1,iint+1))
	    endif
	  else
		if ( si .lt. 0. ) jint = max0(1,min0(njp1,jint-1))
		if ( sj .lt. 0. ) iint = max0(1,min0(nim1,iint-1))
	  endif
	  if ( iold .eq. iint .and. jold .eq. jint ) go to 397
  399	continue
c
	  write (6,*) 'tried maximum number of times but could not find',
     $		'  cell, j,k = ',j,k
	go to 395
  397	continue
	if ( ntries .le. 1 ) then
	   ntries = ntries+1
	   iint = ni
	   jint = jouter
	   go to 417
        else
	     write (6,*) 'not in cell, but no update, i,j = ',iint,jint
     $				,' i1,j1 =  ',j,k
	     write (6,*) 'xzero,yzero',xzero,yzero
	endif
  395	continue
cc
c
 1177   continue
c
cc
c
c	figure out where in the cell we are
c
c
	  ii = eta1(k)
	  jj = chi1(k)

c         bkress (temp. hack)
          if (ii.eq.0.or.jj.eq.0) then
             print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
             print *,'k,ii,jj=',k,ii,jj
             print *,'x,y,z=',points(1,k),points(2,k),points(3,k)
             ii=1
             jj=1
             stop 
             print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
          endif

c        		regular quadrilaterals
	  dx = xzero - x2mid(ii,jj)
	  dy = yzero - y2mid(ii,jj)
	  xeta = x2mid(ii+1,jj)-x2mid(ii,jj)
	  xchi = x2mid(ii,jj+1)-x2mid(ii,jj)
	  xetachi = x2mid(ii+1,jj+1)+x2mid(ii,jj)-x2mid(ii+1,jj)-
     $		     x2mid(ii,jj+1)
	  yeta = y2mid(ii+1,jj)-y2mid(ii,jj)
	  ychi = y2mid(ii,jj+1)-y2mid(ii,jj)
	  yetachi = y2mid(ii+1,jj+1)+y2mid(ii,jj)-y2mid(ii+1,jj)-
     $		     y2mid(ii,jj+1)
	  aa = yeta*xetachi-xeta*yetachi
	  bb = xchi*yeta-xeta*ychi +yetachi*dx - xetachi*dy
	  cc = ychi*dx - xchi*dy
	  disc = sqrt(bb**2-4*aa*cc)/(2.*aa)
	  bt = -bb/(2.*aa)
	  if ( bt + disc .ge. 0. .and. bt + disc .le. 1. ) then
		etainc = bt +disc
	  else
		etainc = bt - disc
	  endif
	  if ( abs(xchi) .ge. abs(ychi) ) then 
	    chiinc = (dx -xeta*etainc)/(xchi+xetachi*etainc)
	  else
	    chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	  endif
	  eta1(k) = float(ii) + etainc
	  chi1(k) = float(jj) + chiinc
	  if ( etainc .gt. 1. .or. etainc .lt. 0. .or.
     $		chiinc .gt. 1. .or. chiinc .lt. 0. ) then
	     write (6,*) 'bad solution for position within cell',
     $		'   i1,j1,i2,j2,eta,chi  ',j,k,ii,jj,etainc,chiinc
		et1 = bt+disc
	  ch1 = (dx -xeta*et1)/(xchi+xetachi*et1)
		et2 = bt-disc
	  ch2 = (dx -xeta*et2)/(xchi+xetachi*et2)
	     write (6,*) 'bt,disc,+ and -',bt,disc,et1,ch1,et2,ch2
	   endif
c
	   i00 = eta1(k)
	   j00 = chi1(k)
	   k00 = xi1(k)
	   deta = eta1(k)-float(i00)
	   dchi = chi1(k)-float(j00)
	   dxi  = xi1(k)-float(k00)
	   icube(k) = (k00-1)*ni*njp2+(j00-1)*ni+i00
	   a(k,1) = (1.-deta)*(1.-dchi)*(1.-dxi)
	   a(k,2) = (1.-deta)*dchi*(1.-dxi)
	   a(k,3) = (1.-deta)*(1.-dchi)*dxi
	   a(k,4) = (1.-deta)*dchi*dxi
	   a(k,5) = deta*(1.-dchi)*(1.-dxi)
	   a(k,6) = deta*dchi*(1.-dxi)
	   a(k,7) = deta*(1.-dchi)*dxi
	   a(k,8) = deta*dchi*dxi
c
c
        ELSE
c
*             ikk = mod(k,30) + 1
*             jkk = mod(k/30,30) + 1
*             kii = mod(k/900,30) + 1
*             write (6,*) 'zeroing a at',k,ikk,jkk,kii
c just to be safe give eta1 and chi1 values to things in the hole
	     eta1(k) = 2.0
	     chi1(k) = 2.0
	   i00 = eta1(k)
	   j00 = chi1(k)
	   k00 = xi1(k)
	   icube(k) = (k00-1)*ni*njp2+(j00-1)*ni+i00
	     do 601 kkkl=1,8
	    a(k,kkkl) = 0.
  601	     continue
c
c
        ENDIF
c
c
  400	    continue
c
c

	return
	end

	function cross(x1,x2)
	dimension x1(2),x2(2)
	real cross !bkress
	cross = x1(1)*x2(2) - x2(1)*x1(2)
	return
	end

	integer function isamin(n,x,incx)
	dimension x(n)
c
	ind = 1
	xind = abs(x(1))
	do 200 i=2,n,incx
	  if ( abs(x(i)) .lt. xind ) then
		ind = i
		xind = abs(x(i))
	  endif
  200	continue
  
c
	isamin = ind
	return
	end
	  







