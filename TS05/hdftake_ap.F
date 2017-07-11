      INTEGER FUNCTION hdftake(first,fname,mstep,xlim,interp,
     $     allcart,points, npts, stime)
c 
c This program reads HDF dump files from the Clock Dial simulation.
c

c USER NOTES:
c This program reads an LFM HDF output file (fname) and puts the mhd fields on 
c a user specified grid (points). The 14 mhd var. are output in allcart(14,i)
c with each i corresponding to (x=points(1,i), y=points(2,i), z=points(3,i)),
c where i = 1 to npts. If first is true, hdftake opens fname and reads in the 
c earth's dipole field data and the mesh metric data (i.e. the first two data 
c sets in the LFM dump files) then reads up to and reads in mstep. If first is 
c false then hdftake assumes that the LFM file is open and reads up to and 
c reads in mstep. xlim appears to do nothing. Set interp=.true. 
c
c	the mhd variables are returned in allcart(14,i) with
c  1 density units particles/cm^3
c	2,3,4 vx,vy,vz km/s
c	5 pressure kev/cm^3
c	6,7,8 Bx,By,Bz nanotesla 
c	9,10,11 Ex,Ey,Ez volts/meter
c	12,13,14 Jx,Jy,Jz guass/cm
c  THESE ARE IN SM COORDINATES
c  - bkress

c df.inc


      integer DFOPEN,DFCLOSE,DFPUT,DFIOPEN
      integer DFERRNO, DFGET
      integer DSSDIMS,DSPDATA,DSADATA
      integer DSGDIMS,DSGDATA

c param.inc
#include "param.inc"
#include "inparam.inc"
c var.inc
#include "var.inc"
#include "pavar.inc"
c meter.inc
      COMMON /METER/
     $     X(NIP1,NJP1,NKP1),Y(NIP1,NJP1,NKP1),Z(NIP1,NJP1,NKP1),
     $     VOLUME(NI,NJ,NK),FACEqI(NIP1,NJ,NK),FACEqJ(NI,NJP1,NK),
     $     X2(NIP1,NJP1),Y2(NIP1,NJP1),PHI(NKP1),
     $     COSPHI(NKP1),SINPHI(NKP1),CPHIqM(NK),SPHIqM(NK),DARC(NK),
     $     FACEqK(NI,NJ,NKP1),COSARC(NK),
     $     EDGEqK(NIP1,NJP1),EDGEqI(NIP1,NJ),EDGEqJ(NI,NJP1),
     $     YBARI(NIP1,NJP1),YBARJ(NIP1,NJP1),DX2I(NIP1,NJ),
     $     DY2I(NIP1,NJ),DX2J(NI,NJP1),DY2J(NI,NJP1),
     $     NRINGS,NAVER(NJ),VOLFACT(NJ,NK)
      parameter (n2dread = 2*nip1*njp1 + nkp1)
      dimension read2d(n2dread)
      equivalence (read2d,x2)

c bzero.inc
	COMMON /BZERO/
     $     BNqFI(NOUTP1,NJ,NK),BsqqI(NOUTp1,NJ,NK),
     $     BZqnI(NOUTP1,NJ,NK),BXqnI(NOUTP1,NJ,NK),BYqnI(NOUTP1,NJ,NK),
     $     BZqFI(NOUTP1,NJ,NK),BXqFI(NOUTP1,NJ,NK),BYqFI(NOUTP1,NJ,NK),
     $     BxqK(NOUTP1,NJP1,NK),ByqK(NOUTP1,NJP1,NK),
     $	      BZQK(NOUTP1,NJP1,NK),
     $     BNqFJ(NOUT,NJP1,NK),BsqqJ(NOUT,NJp1,NK),
     $     BZqnJ(NOUT,NJP1,NK),BXqnJ(NOUT,NJP1,NK),BYqnJ(NOUT,NJP1,NK),
     $     BZqFJ(NOUT,NJP1,NK),BXqFJ(NOUT,NJP1,NK),BYqFJ(NOUT,NJP1,NK),
     $     BxqI(NOUT,NJP1,NKP1),ByqI(NOUT,NJP1,NKP1),
     $	      BZQI(NOUT,NJP1,NKP1),
     $     BNqFK(NOUT,NJ,NKP1),BsqqK(NOUT,NJ,NKp1),
     $     BZqnK(NOUT,NJ,NKP1),BXqnK(NOUT,NJ,NKP1),BYqnK(NOUT,NJ,NKP1),
     $     BZqFK(NOUT,NJ,NKP1),BXqFK(NOUT,NJ,NKP1),BYqFK(NOUT,NJ,NKP1),
     $     BxqJ(NOUTP1,NJ,NKP1),ByqJ(NOUTP1,NJ,NKP1),
     $	      BZQJ(NOUTP1,NJ,NKP1)
       parameter (nbtot = 8*noutp1*nj*nk + 3*noutp1*njp1*nk +
     $      8*nout*njp1*nk + 3*nout*njp1*nkp1 +
     $ 		8*nout*nj*nkp1 + 3*noutp1*nj*nkp1 )
           dimension bzeros(nbtot)
           equivalence (bnqfi,bzeros)

c boundx.inc
	COMMON /WIND/
     $     RHOqIN,VXqIN,VYqIN,VZqIN,CqIN,BXqIN,BYqIN,BZqIN,
     $     RHOq0,Cq0,
     $     MJIN,MJINP1,mjinp2,MJWALL,MJWP1,mjwp2,
     $     DUMMYIN(10,NJ2,NK,9)

#include "iongeo.inc"
c
c	mjin = the number of inflowing j-faces.
c       mjinp1 = the j-value of the front corner?
c	mjwall = the j-value of the rear corner?
c

c
      integer ishape(3),idims(4)
c     added interp here - bkress
      logical first,intfirst,interp
      character*80 fname
      parameter(lentext=80)
      parameter( nctot = nci+ncj+nck)
      character*(lentext) text
      real jd,stime
      real newvar(18,ni,njp2,nkp1)
*     real newcoor(ni,njp2,nkp1,3)
      dimension allcart(15,npts),xlim(6)
      dimension varion(4,njp1,nk2p3),coordion(njp1,nk2p3,2)
      dimension points(3,npts), curtmp(3,ni,njp2,nk)
      dimension zrot(3,3),xrot(3,3),arot(3,3)
      dimension xd(3),vd(3),bd(3),ed(3),jd(3)
      dimension cx(3),cy(3),cz(3),et(3)
c     dimension eoldi(ni,2)

*      print *,'noutp1,nj,nk=',noutp1,nj,nk
c
      hdftake = 0 ! no error
      leng = len(fname)

c     changed the following to if-else loop - bkress
*      if (first) intfirst=.true.
      if (first) then
         intfirst=.true.
      else
         intfirst=.false.
      endif


c*******************************************************************
c*******************************************************************
c*******************************************************************

      if (first) then

         ierr = DFCLOSE(ifile)
*	      print *,'hdf file closed'

c
*      print *, ' Read in the WIND data: mjin,mjwall '
c     set the boundx.inc parameters originally from MESH().
      nshell = nip1 ! nip1 ?
      ntheta = njp1 ! njp1 ?
      nangle = nkp1
      mjin   = 4    ! number of j-faces in the front?
      mjwall = 11 ! back corner ?

c     Read in the Earth's magnetic dipole field data: "bzero.inc" 
*      print *, ' Read in the magnetic dipole field data: "bzero.inc"'
      ierr = DSGDIMS(fname,irank,ishape,3) ! DFSDgetdims()

*         print *,'Read in the magnetic dipole field data'
*         print *,'fname=',fname
*         print *,'irank=',irank
*         print *,'ishape=',ishape

      if( ierr .ne. 0 ) then
         print *,' DFSDgetdims() = ',ierr
         print *,' DFerror = ',DFERRNO()
         hdftake = -1 ! Can not find a data set.
         return
      endif
      if( ishape(1) .ne. nbtot ) then
         print *,' Incompatible nbtot,ishape(1)=',nbtot,ishape(1)

         stop

      endif
      ierr = DSGDATA(fname,irank,ishape,bzeros) ! DFSDgetdata()
      if( ierr .ne. 0 ) then
         print *,' DFSDgetdata(BZEROS) = ',ierr
         print *,' DFerror = ',DFERRNO()
         stop
      endif

c     Write the mesh metric data: "meter.inc" 
*      print *, ' Read in the mesh metric data: "meter.inc"'
      ierr = DSGDIMS(fname,irank,ishape,3) ! DFSDsetdims()

*         print *
*         print *,'fname=',fname
*         print *,'irank=',irank
*         print *,'ishape=',ishape

      if( ierr .ne. 0 ) then
         print *,' DFSDgetdims() = ',ierr
         print *,' DFerror = ',DFERRNO()
         stop
      endif
      if( ishape(1) .ne. n2dread ) then
         print *,' Incompatible n2dread,ishape(1)=',n2dread,ishape(1)

         stop

      endif
      ierr = DSGDATA(fname,irank,ishape,read2d) ! DFSDgetdata()
      if( ierr .ne. 0 ) then
         print *,' DFSDgetdata(READ2D) = ',ierr
         print *,' DFerror = ',DFERRNO()
         stop
      endif

c
      ENDIF
c

c*******************************************************************
c*******************************************************************
c*******************************************************************

c     Process var.inc
      number = 0
 707  continue
*      print *, ' Read in the MHD state data: "var.inc"'
      varall(1) = real(lstep)
      ierr = DSGDIMS(fname,irank,ishape,3) ! DFSDsetdims()

*         print *
*         print *,'fname=',fname
*         print *,'irank=',irank
*         print *,'ishape=',ishape

      if( ierr .ne. 0 ) then
         print *,' DFSDgetdims() = ',ierr
         print *,' DFerror = ',DFERRNO()
         hdftake = -1 ! Can not find a data set.
         return
      endif
      if( ishape(1) .ne. nvarall ) then
         print *,' Incompatible nvarall,ishape(1)=',nvarall,ishape(1)
         if ( ishape(1) .eq. (nvar1+2) ) then
            write (6,*) ' appear just to be missing electric field data'
     $           ,', will continue'
         else

         stop
         endif
      endif
      ierr = DSGDATA(fname,irank,ishape,varall) ! DFSDgetdata()
      if( ierr .ne. 0 ) then
         print *,' DFSDgetdata(VARALL) = ',ierr
         print *,' DFerror = ',DFERRNO()
         stop
      endif
      lstep = int(varall(1))
      print *,'read lstep =',lstep
c

c     get dimensions of pa data set
      ierr = DSGDIMS(fname,irank,ishape,3) ! DFSDsetdims()
*         print *
*         print *,'fname=',fname
*         print *,'irank=',irank
*         print *,'ishape=',ishape
      if( ierr .ne. 0 ) then
         print *,' DFSDgetdims() = ',ierr
         print *,' DFerror = ',DFERRNO()
         hdftake = -1 ! Can not find a data set.
         return
      endif
      if( ishape(1) .ne. ntot ) then
         print *,' Incompatible ntot,ishape(1)=',ntot,ishape(1)
         stop
      endif
c     get pa data     
      ierr = DSGDATA(fname,irank,ishape,pavar) ! DFSDgetdata()
      if( ierr .ne. 0 ) then
         print *,' DFSDgetdata(pavar) = ',ierr
         print *,' DFerror = ',DFERRNO()
         stop
      endif



      number = number+1
      if ( mstep .gt. lstep .and. ierr .eq. 0 .and. number .lt. 500)
     $     go to 707

        if (first) then
	call metric
        endif


c*******************************************************************
c*******************************************************************
c*******************************************************************

      stime = varall(2)
      print *,'In hdftake: read lfm step = ',lstep
      print *,'time(sec) = ',stime

c
c
c
c
c    start writing the new data
c       1 - density units particles/cm^3
c	2,3,4 vx,vy,vz units kms
c	5 pressure kev/cm^3
c	6,7,8 Bx,By,Bz units of nanotesla
c	9,10,11 x,y,z units of Re
c
	do 700 k=1,nk
	do 700 j=1,nj
	do 700 i=1,ni
	  newvar(1,i,j+1,k) = rho(i,j,k)/2.1e-24
	  newvar(2,i,j+1,k) = vx(i,j,k)/1.e5
	  newvar(3,i,j+1,k) = vy(i,j,k)/1.e5
	  newvar(4,i,j+1,k) = vz(i,j,k)/1.e5
	  newvar(5,i,j+1,k) = 0.6*rho(i,j,k)*c(i,j,k)**2/1.6e-9
	  newvar(6,i,j+1,k) = (bx(i,j,k)+
     $		0.5*(bxqnk(i,j,k)+bxqnk(i,j,k+1)))*1.e5
	  newvar(7,i,j+1,k) = (by(i,j,k)+
     $		0.5*(byqnk(i,j,k)+byqnk(i,j,k+1)))*1.e5
	  newvar(8,i,j+1,k) = (bz(i,j,k)+
     $		0.5*(bzqnk(i,j,k)+bzqnk(i,j,k+1)))*1.e5

	 newvar(15,i,j+1,k) =pa(i,j,k)

	 newvar(16,i,j+1,k) =
     $		 0.125*(x(i,j,k)+x(i+1,j,k)+x(i,j+1,k)+x(i+1,j+1,k)+
     $	  x(i,j,k+1)+x(i,j+1,k+1)+x(i+1,j,k+1)+x(i+1,j+1,k+1))/rearth
	 newvar(17,i,j+1,k) = 
     $		0.125*(y(i,j,k)+y(i+1,j,k)+y(i,j+1,k)+y(i+1,j+1,k)+
     $	  y(i,j,k+1)+y(i,j+1,k+1)+y(i+1,j,k+1)+y(i+1,j+1,k+1))/rearth
	 newvar(18,i,j+1,k) = 
     $		0.125*(z(i,j,k)+z(i+1,j,k)+z(i,j+1,k)+z(i+1,j+1,k)+
     $	  z(i,j,k+1)+z(i,j+1,k+1)+z(i+1,j,k+1)+z(i+1,j+1,k+1))/rearth
  700	continue
c
c  get cartesian E's at regular grid points
c        do 704 i=1,ni
c           enew1 = (eoldi(i,1)-etoti(i,1,1))*0.01
c           enew2 = (eoldi(i,2)-etoti(i,nj,1))*0.01
c           do 704 k=1,nkp1
c              etoti(i,1,k) = enew1
c              etoti(i,nj,k) = enew2
c 704           continue
        do 710 k=1,nk
           do 710 j=1,nj
              do 710 i=1,ni
                 e_i = 0.25*(etoti(i,j,k)+etoti(i,j,k+1) +
     $                 etoti(i,j+1,k)+etoti(i,j+1,k+1))
                 e_j = 0.25*(etotj(i,j,k)+etotj(i,j,k+1) +
     $                 etotj(i+1,j,k)+etotj(i+1,j,k+1))
                 e_k = 0.25*(etotk(i,j,k)+
     $                etotk(i+1,j,k) +
     $                etotk(i,j+1,k)+
     $                etotk(i+1,j+1,k))
                 x_i = 0.25*(
     $           (x(i+1,j,k)+x(i+1,j+1,k)+x(i+1,j,k+1)+x(i+1,j+1,k+1))
     $        -  (x(i,j,k)+x(i,j+1,k)+x(i,j,k+1)+x(i,j+1,k+1)))
                 y_i = 0.25*(
     $           (y(i+1,j,k)+y(i+1,j+1,k)+y(i+1,j,k+1)+y(i+1,j+1,k+1))
     $        -  (y(i,j,k)+y(i,j+1,k)+y(i,j,k+1)+y(i,j+1,k+1)))
                 z_i = 0.25*(
     $           (z(i+1,j,k)+z(i+1,j+1,k)+z(i+1,j,k+1)+z(i+1,j+1,k+1))
     $        -  (z(i,j,k)+z(i,j+1,k)+z(i,j,k+1)+z(i,j+1,k+1)))
                 x_j = 0.25*(
     $           (x(i,j+1,k)+x(i+1,j+1,k)+x(i,j+1,k+1)+x(i+1,j+1,k+1))
     $        -  (x(i,j,k)+x(i+1,j,k)+x(i,j,k+1)+x(i+1,j,k+1)))
                 y_j = 0.25*(
     $           (y(i,j+1,k)+y(i+1,j+1,k)+y(i,j+1,k+1)+y(i+1,j+1,k+1))
     $        -  (y(i,j,k)+y(i+1,j,k)+y(i,j,k+1)+y(i+1,j,k+1)))
                 z_j = 0.25*(
     $           (z(i,j+1,k)+z(i+1,j+1,k)+z(i,j+1,k+1)+z(i+1,j+1,k+1))
     $        -  (z(i,j,k)+z(i+1,j,k)+z(i,j,k+1)+z(i+1,j,k+1)))
                 x_k = 0.25*(
     $           (x(i,j,k+1)+x(i+1,j,k+1)+x(i,j+1,k+1)+x(i+1,j+1,k+1))
     $        -  (x(i,j,k)+x(i,j+1,k)+x(i+1,j,k)+x(i+1,j+1,k)))
                 y_k = 0.25*(
     $           (y(i,j,k+1)+y(i+1,j,k+1)+y(i,j+1,k+1)+y(i+1,j+1,k+1))
     $        -  (y(i,j,k)+y(i,j+1,k)+y(i+1,j,k)+y(i+1,j+1,k)))
                 z_k = 0.25*(
     $           (z(i,j,k+1)+z(i+1,j,k+1)+z(i,j+1,k+1)+z(i+1,j+1,k+1))
     $        -  (z(i,j,k)+z(i,j+1,k)+z(i+1,j,k)+z(i+1,j+1,k)))
                 cx(1) = x_i
                 cx(2) = x_j
                 cx(3) = x_k
                 cy(1) = y_i
                 cy(2) = y_j
                 cy(3) = y_k
                 cz(1) = z_i
                 cz(2) = z_j
                 cz(3) = z_k
                 et(1) = e_i
                 et(2) = e_j
                 et(3) = e_k
                 det = 1.e-6/triple(cx,cy,cz)
                 newvar(9,i,j+1,k) = triple(et,cy,cz)*det
                 newvar(10,i,j+1,k) = triple(cx,et,cz)*det
                 newvar(11,i,j+1,k) = triple(cx,cy,et)*det
 710    continue
c
c
c        get  J's at cell centers
c
        call flux2ave
        do k=1,nk
        do j=1,nj
        do i=1,nip1
           bi(i,j,k) = bi(i,j,k)+bnqfi(i,j,k)
        enddo
        enddo
        enddo
c
        do k=1,nk
        do j=1,njp1
        do i=1,ni
           bj(i,j,k) = bj(i,j,k)+bnqfj(i,j,k)
        enddo
        enddo
        enddo
c
        do k=1,nkp1
        do j=1,nj
        do i=1,ni
           bk(i,j,k) = bk(i,j,k)+bnqfk(i,j,k)
        enddo
        enddo
        enddo

c       currents, i.e. 12,13,14
        call bzz(bx,by,bz,bi,bj,bk)
        call jcalc2(curtmp,bx,by,bz)
        do l=1,3
        do k=1,nk
           do j=1,njp2
              do i=1,ni
                 newvar(11+l,i,j,k) = curtmp(l,i,j,k)
                 enddo
           enddo
        enddo
        enddo



c
c      for all MHD variables
***	do 810 l=1,17
      do 810 l=1,18
	do 810 j=1,njp2,njp1
	do 810 i=1,ni
	  var_av=0.0
	  jj0 = max(2,min(njp1,j))
	do 800 k=1,nk
	  var_av = var_av + newvar(l,i,jj0,k)
  800	continue
	var_av = var_av/float(nk)
	do 805 k=1,nk
	  newvar(l,i,j,k) = var_av
  805	continue
  810	continue
c
c
	do 840 j=1,njp2
	do 840 i=1,ni

c      for all MHD variables
***	do 840 l=1,17
      do 840 l=1,18
	  newvar(l,i,j,nkp1) = newvar(l,i,j,1)
  840	continue
c

c      for XYZ coors.
*	do 860 l=15,17
*	do 860 k=1,nkp1
*	do 860 j=1,njp2
*	do 860 i=1,ni
*	  newcoor(i,j,k,l-14) = newvar(l,i,j,k)
*  860	continue
c         
c
        if (interp) then
           call cart1d(newvar,points,15,18,npts,allcart,intfirst)
        endif
c

      RETURN
      END

c
        function triple(x,y,z)
        dimension x(3),y(3),z(3)
        dimension dum(3)
        call cross3(y,z,dum)
        triple = dot(x,dum)
        return
        end
c
c
        function dot(x,y)
        dimension x(3),y(3)
        dot = x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
        return
        end
c
c
        subroutine cross3(x,y,z)
        dimension x(3),y(3),z(3)
        z(1) = x(2)*y(3) - x(3)*y(2)
        z(2) = x(3)*y(1) - x(1)*y(3)
        z(3) = X(1)*y(2) - x(2)*y(1)
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

	SUBROUTINE INT2STR(CAOUT,IVALUE,IBEGIN,IEND)

	character*(*) caout
	INTEGER iValue, iBegin, iEnd
	INTEGER i,iStart,iLast
	character*20 caDum


	WRITE(caDum,10) iValue
 10     FORMAT(i8)

	iStart = 1
	DO while( caDum(iStart:iStart) .eq. ' ' )
	  iStart = iStart + 1
	ENDDO
*#ifdef AIX
*	iLast = LEN_TRIM(caDum)
*#endif
*#ifndef AIX
 	iLast = LNBLNK(caDum)
*#endif

	caOut(1:iLast-iStart+1) = caDum(iStart:iLast)
        iBegin = 1
	iEnd = iLast - iStart + 1

	RETURN
	END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	SUBROUTINE FLT2STR(CAOUT,RVALUE,IBEGIN,IEND)

	character*(*) caout
	real rValue
	INTEGER iBegin, iEnd
	INTEGER i,iStart,iLast
	character*20 caDum


	WRITE(caDum,10) rValue
 10     format(1pe16.7)

	iStart = 1
	DO while( caDum(iStart:iStart) .eq. ' ' )
	  iStart = iStart + 1
	ENDDO
#ifdef AIX
	iLast = LEN_TRIM(caDum)
#endif
#ifndef AIX
 	  iLast = LNBLNK(caDum)
#endif

	
	caOut(1:iLast-iStart+1) = caDum(iStart:iLast)
        iBegin = 1
	iEnd = iLast - iStart + 1

	RETURN
	END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	SUBROUTINE RSTSTR(CAOUT)
	
	character*(*) caOut

	iLen=len(caOut)
	DO i=1,iLen
#ifdef AIX
  	  caOut(i:i)='\0'
#endif
#if defined (SGI) || defined (LINUX)
          caOut(i:i)=achar(0)
#endif
	ENDDO

	RETURN
	END
