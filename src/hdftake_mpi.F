
	INTEGER FUNCTION hdftake(fname,xlim,interp,allcart,points,npts,
     &  stime)

	CHARACTER*(*) fname
        DIMENSION xlim(6)
        LOGICAL interp
        DIMENSION allcart(14,npts),points(3,npts)
        INTEGER npts
        real stime

#	include "hdf.inc"
#	include "help.inc"
#	include "param.inc"
#	include "meter.inc"
#	include "var.inc"
#	include "ionC90.inc"
#	include "bzero.inc"
#	include "boundx.inc"

	INTEGER*4 iFileId, iSecId,iSecIndex, iaDimId, iAttrId, iRank
	INTEGER iNtype, iNAttr
	INTEGER sfstart, sfend, sfscatt, sfcreate, sfwcdata
	INTEGER sfwdata, sfendacc,sffattr,sfrnatt,sfrcatt
	INTEGER sfn2index,sfginfo,sfrdata,sfselect
 	INTEGER sfsattr, sfrattr
	INTEGER iLen, iStatus 
	INTEGER iFileStart, iFileStop, iFileInterval
	INTEGER iaDim4d(4), iaStart4d(4), iaEdge4d(4),iaStride4d(4)
	INTEGER iaDim3d(3), iaStart3d(3), iaEdge3d(3),iaStride3d(3)
	INTEGER iaDim2d(2), iaStart2d(2), iaEdge2d(2),iaStride2d(2)
	INTEGER iaDim1d(1), iaStart1d(1), iaEdge1d(1),iaStride1d(1)
	CHARACTER*80 caCont, caDim, caSec, caAttr
	CHARACTER*20 caTmp
	LOGICAL InFile,intfirst

        REAL newvar(17,ni,njp2,nkp1),newcoor(ni,njp2,nkp1,3)
        REAL curtmp(3,ni,njp2,nk)
        REAL cx(3),cy(3),cz(3),et(3)

        REAL raVar(nip1,njp1,nkp1)
        CHARACTER caVars(14)*5

c        for no low pass
*        DATA caVars/'rho_','vx_','vy_','vz_','c_','bx_','by_','bz_',
*     $  'bi','bj','bk','ei','ej','ek'/

c       for no low pass
*        DATA caVars/'rho_','vx_','vy_','vz_','c_','bx_','by_','bz_',
*     $  'bi_','bj_','bk_','ei_','ej_','ek_'/

c       (alt. var. naming convention)
c       for low pass
c       NOTE THAT HERE WE HAVE PERTURBATION B ONLY!!
*        DATA caVars/'rho_','vx_','vy_','vz_','c_','avgBx_','avgBy_',
*     $  'avgBz_','bi','bj','bk','avgEi','avgEj','avgEk'/

c       for low pass
c       NOTE THAT HERE WE HAVE PERTURBATION B ONLY!!
        DATA caVars/'rho_','vx_','vy_','vz_','c_','avgBx','avgBy',
     $  'avgBz','bi_','bj_','bk_','avgEi','avgEj','avgEk'/

*        print *,'caVars=',caVars

	iFileId = sfstart(fname,DFACC_RDONLY)
	IF (iFileId .lt. 0 ) THEN
	   WRITE(6,*) 'Unable to open ',fname,' for read!'
	   STOP
	ENDIF
	
  	!Read in Grid info
  	iSecIndex = sfn2index(iFileId,'X_grid')
  	iSecId = sfselect(iFileId, iSecIndex)
  	iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
        write(6,*) iaDim3D
        IF (iStatus .EQ. -1) THEN
          WRITE(6,*) ' HDFTAKE: Failure reading grid info'
          STOP
        ENDIF
  	IF ( nip1 .ne. iaDim3d(1) .or. njp1 .ne. iaDim3d(2)) THEN
  	    WRITE(6,*) ' HDFTAKE: File and code size incompatible'

            print *,'nip1=',nip1
            print *,'iaDim2d(1)=',iaDim2d(1)

  	    STOP
  	ENDIF
  	DO i=1,iRank
  	    iaStride3d(i) = 1
  	    iaStart3d(i) = 0
    	    iaEdge3d(i) = iaDim3d(i)
  	ENDDO
  	iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,X)
        IF (iStatus .EQ. -1)  THEN
            WRITE(6,*) ' HDFTAKE: Failure reading X_grid'
  	    STOP
        ENDIF
  	iStatus = sfendacc(iSecId)
        IF (iStatus .EQ. -1)  THEN
            WRITE(6,*) ' HDFTAKE: Failure closing X_grid'
     	    STOP
  	ENDIF

  	iSecIndex = sfn2index(iFileId,'Y_grid')
  	iSecId = sfselect(iFileId, iSecIndex)
  	iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
        IF (iStatus .EQ. -1) THEN
          WRITE(6,*) ' HDFTAKE: Failure reading grid info'
          STOP
        ENDIF
  	IF ( nip1 .ne. iaDim3d(1) .or. njp1 .ne. iaDim3d(2)) THEN
  	    WRITE(6,*) ' HDFTAKE: File and code size incompatible'

            print *,'nip1=',nip1
            print *,'iaDim2d(1)=',iaDim2d(1)

  	    STOP
  	ENDIF
  	DO i=1,iRank
  	    iaStride3d(i) = 1
  	    iaStart3d(i) = 0
    	    iaEdge3d(i) = iaDim3d(i)
  	ENDDO
  	iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,Y)
        IF (iStatus .EQ. -1)  THEN
            WRITE(6,*) ' HDFTAKE: Failure reading Y_grid'
  	    STOP
        ENDIF
  	iStatus = sfendacc(iSecId)
        IF (iStatus .EQ. -1)  THEN
            WRITE(6,*) ' HDFTAKE: Failure closing Y_grid'
     	    STOP
  	ENDIF

  	iSecIndex = sfn2index(iFileId,'Z_grid')
  	iSecId = sfselect(iFileId, iSecIndex)
  	iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
        IF (iStatus .EQ. -1) THEN
          WRITE(6,*) ' HDFTAKE: Failure reading grid info'
          STOP
        ENDIF
  	IF ( nip1 .ne. iaDim3d(1) .or. njp1 .ne. iaDim3d(2)) THEN
  	    WRITE(6,*) ' HDFTAKE: File and code size incompatible'

            print *,'nip1=',nip1
            print *,'iaDim3d(1)=',iaDim3d(1)

  	    STOP
  	ENDIF
  	DO i=1,iRank
  	    iaStride3d(i) = 1
  	    iaStart3d(i) = 0
    	    iaEdge3d(i) = iaDim3d(i)
  	ENDDO
  	iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,Z)
        IF (iStatus .EQ. -1)  THEN
            WRITE(6,*) ' HDFTAKE: Failure reading Z_grid'
  	    STOP
        ENDIF
  	iStatus = sfendacc(iSecId)
        IF (iStatus .EQ. -1)  THEN
            WRITE(6,*) ' HDFTAKE: Failure closing Z_grid'
     	    STOP
  	ENDIF
        DO k = 1,nkp1
           phi(k) = acos(y(1,1,k)/y(1,1,1))
           phi(k) = atan2(z(1,1,k),y(1,1,k))
           if( phi(k) .lt. 0.0) then
             phi(k) = phi(k)+2.0*3.141592653589793
           endif
           cosphi(k) = cos(phi(k))
           sinphi(k) = sin(phi(k))
        ENDDO           
        DO k = 1,nkp1
          DO j = 1,njp1
            DO i = 1,nip1
               x2(i,j) = X(i,j,1)
               y2(i,j) = Y(i,j,k)/(cosphi(k))
            ENDDO
          ENDDO
        ENDDO

c       bkress - needed hack
        phi(nkp1)=2.0*3.141592653589793

*********************************************************
*c       bkress
*        DO k = 1,nkp1
*           print *,'phi=',phi(k)
*        ENDDO 
*c       at k=0, phi=0
*c       at k=nkp1, phi=1.74845553E-07
*
*      OPEN(UNIT=82,FILE='GRIDOUT',FORM='FORMATTED',STATUS='OLD')
*      READ(82,102) nshell,ntheta,mjin,mjwall
*  102 FORMAT(4i6)
*
*      print *
*      print *,'nshell,ntheta,mjin,mjwall=',nshell,ntheta,mjin,mjwall
*
*      IF ( nshell .EQ. nip1 .AND. nthetA .EQ. njp1) GO TO 77
*      WRITE (6,151)
*      WRITE (9,151)
*  151 FORMAT(1X,'NOT A PROPER DATA FILE FOR THE COORDINATES')
*      STOP
*   77 CONTINUE
*c
*      READ(82,101) (phi(I),I=1,nkp1)
*      READ(82,101) ((x2(I,J),I=1,nip1),J=1,njp1)
*      READ(82,101) ((y2(I,J),I=1,nip1),J=1,njp1)
*  101 FORMAT(5e16.7)
*      CLOSE(82)
*
**        print *
*        DO k = 1,nkp1
*           phi(k)=phi(k)*3.141592653589793/180.
*           print *,'phi=',phi(k)
*        ENDDO 
*c       at k=0, phi=0
*c       at k=nkp1, phi=6.2831855 (otherwise same as above)
*
*        DO j = 1,njp1
*          DO i = 1,nip1
*             x2(i,j) = x2(i,j)*rearth
*             y2(i,j) = y2(i,j)*rearth
*          ENDDO
*        ENDDO
*
c       bkress
*           print *
*           do i=1,nip1
*             print *,x(i,6,8)/6378.137e5
*             print *,y(1,1,1)
*           enddo
*           stop
*
*********************************************************

          !Process some grid related quantities
c
*          WRITE (6,*) ' HDFTAKE: mjin and mjwall',mjin,mjwall
          dxm = x2(nip1,2) - x2(nip1,1)
          dym = y2(nip1,2) - y2(nip1,1)
          dlen = sqrt(dxm**2 + dym**2)
          dxm = dxm/dlen
          dym = dym/dlen
          ncorns = 1
          DO 22 j=2,nj
            dxp = x2(nip1,j+1) - x2(nip1,j)
            dyp = y2(nip1,j+1) - y2(nip1,j)
            dlen = sqrt(dxp**2 + dyp**2)
            dxp = dxp/dlen
            dyp = dyp/dlen
            IF ( dxm*dyp-dxp*dym .GT. 0.8 ) THEN
              IF ( ncorns .EQ. 1 ) THEN
                mjin = j-1
                ncorns =2
              ENDIf
            ELSEIF ( dxm*dxp + dym*dyp .LT. -0.9 ) THEN
              IF (ncorns .EQ. 2 ) THEN
                ncorns = 3
                mjwall = j-1
              ENDIF
            ENDIF
   22     CONTINUE
*          WRITE (6,*) ' HDFTAKE: mjin and mjwall',mjin,mjwall

          ! read time
  	  iAttrId = sffattr(iFileId,'time')
c	  iStatus = afgainfo(iFileId,iAttrId,'time',data_type,n_values)
  	  iStatus = sfrnatt(iFileId,iAttrId,time)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading time '
  	    STOP
          ENDIF
          stime = time

  	  ! Read in Data for timestep
          DO l = 1,14
*            print *,'reading #,caVars(#)',l,caVars(l)
    	    iSecIndex = sfn2index(iFileId,caVars(l))
  	    iSecId = sfselect(iFileId, iSecIndex)
  	    iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
            IF (iStatus .EQ. -1) THEN
  	      WRITE(6,*) ' HDFTKAE: Failure geting step info'
  	      STOP
            ENDIF
  	    DO i=1,iRank
  	      iaStride3d(i) = 1
  	      iaStart3d(i) = 0
  	      iaEdge3d(i) = iaDim3d(i)
  	    ENDDO
  	    iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,raVar)
            IF (iStatus .EQ. -1) THEN
  	       WRITE(6,*) 'Failure reading step'
  	       STOP
            ENDIF
  	    iStatus = sfendacc(iSecId)
            IF (iStatus .EQ. -1) THEN
  	      WRITE(6,*) ' HDFTAKE: Failure closing step'
  	      STOP
            ENDIF
            if (l .eq. 1) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,ni
                 rho(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
           if (l .eq. 2) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,ni
                 vx(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 3) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,ni
                 vy(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 4) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,ni
                 vz(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 5) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,ni
                 c(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 6) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,ni
                 bx(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 7) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,ni
                 by(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 8) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,ni
                 bz(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 9) THEN
              DO k = 1,nk
               DO j = 1,nj
                DO i = 1,nip1
                 bi(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 10) THEN
              DO k = 1,nk
               DO j = 1,njp1
                DO i = 1,ni
                 bj(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 11) THEN
              DO k = 1,nkp1
               DO j = 1,nj
                DO i = 1,ni
                 bk(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 12) THEN
              DO k = 1,nkp1
               DO j = 1,njp1
                DO i = 1,ni
                 etoti(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 13) THEN
              DO k = 1,nkp1
               DO j = 1,nj
                DO i = 1,nip1
                 etotj(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
            if (l .eq. 14) THEN
              DO k = 1,nk
               DO j = 1,njp1
                DO i = 1,nip1
                 etotk(i,j,k) = raVar(i,j,k)
                ENDDO
               ENDDO
              ENDDO
           ENDIF
          ENDDO
        HDFTAKE = 1
	iStatus = sfend(iFileId)
        if (iStatus .EQ. -1) then
	  WRITE(6,*) ' HDFTAKE: Failure closing file'
	  STOP
        ENDIF
* 	WRITE(6,*) ' HDFTAKE: Done.'

c
c Now we do some processing to create cartiesan data in the correct units
c

        CALL metric
        
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
	  newvar(6,i,j+1,k) = bx(i,j,k)*1.e5
	  newvar(7,i,j+1,k) = by(i,j,k)*1.e5
	  newvar(8,i,j+1,k) = bz(i,j,k)*1.e5
	 newvar(15,i,j+1,k) =
     $		 0.125*(x(i,j,k)+x(i+1,j,k)+x(i,j+1,k)+x(i+1,j+1,k)+
     $	  x(i,j,k+1)+x(i,j+1,k+1)+x(i+1,j,k+1)+x(i+1,j+1,k+1))/rearth
	 newvar(16,i,j+1,k) = 
     $		0.125*(y(i,j,k)+y(i+1,j,k)+y(i,j+1,k)+y(i+1,j+1,k)+
     $	  y(i,j,k+1)+y(i,j+1,k+1)+y(i+1,j,k+1)+y(i+1,j+1,k+1))/rearth
	 newvar(17,i,j+1,k) = 
     $		0.125*(z(i,j,k)+z(i+1,j,k)+z(i,j+1,k)+z(i+1,j+1,k)+
     $	  z(i,j,k+1)+z(i,j+1,k+1)+z(i+1,j,k+1)+z(i+1,j+1,k+1))/rearth
  700	continue
c
c  get cartesian Es at regular grid points
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
c        get  Js at cell centers
c
        call flux2ave
        do k=1,nk
        do j=1,nj
        do i=1,nip1
           bi(i,j,k) = bi(i,j,k)
        enddo
        enddo
        enddo
c
        do k=1,nk
        do j=1,njp1
        do i=1,ni
           bj(i,j,k) = bj(i,j,k)
        enddo
        enddo
        enddo
c
        do k=1,nkp1
        do j=1,nj
        do i=1,ni
           bk(i,j,k) = bk(i,j,k)
        enddo
        enddo
        enddo
c
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
c
	do 810 l=1,17
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
	do 840 l=1,17
	  newvar(l,i,j,nkp1) = newvar(l,i,j,1)
  840	continue
c
	do 860 l=15,17
	do 860 k=1,nkp1
	do 860 j=1,njp2
	do 860 i=1,ni
	  newcoor(i,j,k,l-14) = newvar(l,i,j,k)
  860	continue

      print *,'i,j,k=',ni,1,1
      print *,'writing field values at x,y,z=',
     &X(ni,1,1)/rearth,Y(ni,1,1)/rearth,Z(ni,1,1)/rearth
      print *,'N(#/cc) is newvar(1,53,1,1)=',newvar(1,53,1,1)
      print *,'vx(km/s) is newvar(2,53,1,1)=',newvar(2,53,1,1)
      print *,'bz(nT) is newvar(8,53,1,1)=',newvar(8,53,1,1)
c
c
        if (interp) then
           intfirst = .true.
           call cart1d(newvar,points,14,17,npts,allcart,intfirst)
        endif

        hdftake = 0
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

*c
*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*c
*
*	SUBROUTINE INT2STR(CAOUT,IVALUE,IBEGIN,IEND)
*
*	character*(*) caout
*	INTEGER iValue, iBegin, iEnd
*	INTEGER i,iStart,iLast
*	character*20 caDum
*
*
*	WRITE(caDum,10) iValue
* 10     FORMAT(i8)
*
*	iStart = 1
*	DO while( caDum(iStart:iStart) .eq. ' ' )
*	  iStart = iStart + 1
*	ENDDO
**#ifdef AIX
**	iLast = LEN_TRIM(caDum)
**#endif
**#ifndef AIX
* 	iLast = LNBLNK(caDum)
**#endif
*
*	caOut(1:iLast-iStart+1) = caDum(iStart:iLast)
*        iBegin = 1
*	iEnd = iLast - iStart + 1
*
*	RETURN
*	END
*c
*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*c
*	SUBROUTINE FLT2STR(CAOUT,RVALUE,IBEGIN,IEND)
*
*	character*(*) caout
*	real rValue
*	INTEGER iBegin, iEnd
*	INTEGER i,iStart,iLast
*	character*20 caDum
*
*
*	WRITE(caDum,10) rValue
* 10     format(1pe16.7)
*
*	iStart = 1
*	DO while( caDum(iStart:iStart) .eq. ' ' )
*	  iStart = iStart + 1
*	ENDDO
*#ifdef AIX
*	iLast = LEN_TRIM(caDum)
*#endif
*#ifndef AIX
* 	  iLast = LNBLNK(caDum)
*#endif
*
*	
*	caOut(1:iLast-iStart+1) = caDum(iStart:iLast)
*        iBegin = 1
*	iEnd = iLast - iStart + 1
*
*	RETURN
*	END
*c
*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*c
*	SUBROUTINE RSTSTR(CAOUT)
*	
*	character*(*) caOut
*
*	iLen=len(caOut)
*	DO i=1,iLen
*#ifdef AIX
*  	  caOut(i:i)='\0'
*#endif
*#if defined (SGI) || defined (LINUX)
*          caOut(i:i)=achar(0)
*#endif
*	ENDDO
*
*	RETURN
*	END
