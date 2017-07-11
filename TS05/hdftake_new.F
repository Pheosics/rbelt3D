	INTEGER FUNCTION hdftake(first,fname,mstep,xlim,interp,
     $                           allcart,points,npts,stime)

	LOGICAL first
	CHARACTER*(*) fname
        INTEGER mstep
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

c
c       moved equlivalence statments to here to prevent OMP2.0 problems
c
c       EQUIVALENCE (read2d,x2)
c       EQUIVALENCE (bnqfi,azeros)

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

        intfirst = .true.
	iFileId = sfstart(fname,DFACC_RDONLY)
	IF (iFileId .lt. 0 ) THEN
	   WRITE(6,*) 'Unable to open ',fname,' for read!'
	   STOP
	ENDIF
	
	!determine who wrote the file

	iSecId = sffattr(iFileId,'written_by')
	iStatus = sfrcatt(iFileId,iSecId,caAttr)
        IF (iStatus .EQ. -1) THEN
           WRITE(6,*) ' HDFTAKE: Failed to find written_by'
           STOP
        ENDIF
	IF (caAttr(1:1) .EQ. 'C') THEN
	   WRITE(6,*) ' HDFTAKE: Unable to restart from C files'
	   STOP
	ELSE
	  ! check to see if the starting step is valid
	  iSecId = sffattr(iFileId,'initial_step')
	  iStatus = sfrnatt(iFileId,iSecId,iFileStart)
          IF (iStatus .EQ. -1) THEN
            WRITE(6,*) ' HDFTAKE: Failed to determine initial_step'
            STOP
          ENDIF
  	  iSecId = sffattr(iFileId,'final_step')
  	  iStatus = sfrnatt(iFileId,iSecId,iFileStop)
          IF (iStatus .EQ. -1) THEN
            WRITE(6,*) ' HDFTAKE: Failed to determine initial_step'
            STOP
          ENDIF
  	  iSecId = sffattr(iFileId,'dump_interval')
  	  iStatus = sfrnatt(iFileId,iSecId,iDumpInterval)
          IF (iStatus .EQ. -1) THEN
            WRITE(6,*) ' HDFTAKE: Failed to determine dump_interval'
            STOP
          ENDIF
#ifdef DEBUG
  	  WRITE(6,*) ' HDFTAKE: IFileStart: ',iFileStart
  	  WRITE(6,*) ' HDFTAKE: IFileStop: ',iFileStop
  	  WRITE(6,*) ' HDFTAKE: IFileDump: ',iDumpInterval
#endif
  	  IF ( mstep .LT. iFileStart ) THEN
  	    InFile = .False.
  	    iGetStep = iFileStart
      	    WRITE(6,*) 
     $     ' HDFTAKE Warning: starting step less than',
     $     ' first step in file' 
  	     WRITE(6,*) ' HDFTAKE: Will return first step in file'
  	  ENDIF
  	  IF ( mstep .GT. iFileStop) THEN
  	    InFile = .False.
  	    iGetStep = iFileStop
  	    WRITE(6,*) 
     $      ' HDFTAKE Warning: starting step greater than',
     $      ' last step in file' 
  	    WRITE(6,*) ' HDFTAKE: Will return last step in file'
          ENDIF
  	  IF ( mstep .GE. iFileStart .AND. mstep .LE. iFileStop) THEN
  	    IF ( MOD(mstep,iDumpInterval) .EQ. 0) THEN
  	      inFile = .TRUE.
  	      iGetStep = mstep
  	    ELSE
  	      InFile = .FALSE.
  	      iGetStep = mstep - MOD(mstep,iDumpInterval)
              WRITE(6,*) ' HDFTAKE Warning: starting step not',
     $                 ' valid dump'
  	      WRITE(6,*) ' HDFTAKE: Will return cloest previous step'
  	    ENDIF
  	  ENDIF
  
  	  WRITE(6,*) ' HDFTAKE: getting step',iGetStep
  
  	  !Read in Grid info
  	  iSecIndex = sfn2index(iFileId,'grid_x2')
  	  iSecId = sfselect(iFileId, iSecIndex)
  	  iStatus = sfginfo(iSecId,caSec,iRank,iaDim2d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
            WRITE(6,*) ' HDFTAKE: Failure reading grid info'
            STOP
          ENDIF
  	  IF ( nip1 .ne. iaDim2d(1) .or. njp1 .ne. iaDim2d(2)) THEN
  	    WRITE(6,*) ' HDFTAKE: File and code size incompatible'

            print *,'nip1=',nip1
            print *,'iaDim2d(1)=',iaDim2d(1)
            print *,'iaDim2d(2)=',iaDim2d(2)

  	    STOP
  	  ENDIF
  	  DO i=1,iRank
  	    iaStride2d(i) = 1
  	    iaStart2d(i) = 0
    	    iaEdge2d(i) = iaDim2d(i)
  	  ENDDO
  	  iStatus = sfrdata(iSecId,iaStart2d,iaStride2d,iaEdge2d,x2)
          IF (iStatus .EQ. -1)  THEN
            WRITE(6,*) ' HDFTAKE: Failure reading x2'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1)  THEN
            WRITE(6,*) ' HDFTAKE: Failure closing x2'
     	    STOP
  	  ENDIF
          iSecIndex = sfn2index(iFileId,'grid_y2')
          iSecId = sfselect(iFileId, iSecIndex)
          iStatus = sfrdata(iSecId,iaStart2d,iaStride2d,iaEdge2d,y2)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading y2'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing y2'
  	    STOP
          ENDIF
#ifdef DEBUG
  	  WRITE(6,*) ' HDFTAKE: grid from inside HDFTAKE'
          DO j=1,5
  	    DO i=1,5
  	     WRITE(6,90) i,j,x2(i,j),y2(i,j)
            ENDDO
          ENDDO
   90     FORMAT(1x,'i: ',i,' j: ',i,' x2:',1e16.7,' y2: ',1e16.7)
#endif
  	  iSecIndex = sfn2index(iFileId,'grid_phi')
  	  iSecId = sfselect(iFileId, iSecIndex)
  	  iStatus = sfginfo(iSecId,caSec,iRank,iaDim1d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE:Failure reading phi info'
  	    STOP
          ENDIF
  	  IF (nkp1 .NE. iaDim1d(1)) THEN
  	    WRITE(6,*) ' HDFTAKE: File and code size incompatible'
  	    STOP
  	  ENDIF
  	  DO i=1,iRank
  	    iaStart1d(i) = 0
  	    iaEdge1d(i) = iaDim1d(i)
  	  ENDDO
  	  iStatus = sfrdata(iSecId,iaStart1d,1,iaEdge1d,phi)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading phi'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing phi'
  	    STOP
          ENDIF
  
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
  
  	  ! Read in Bzero information
  	  iSecIndex = sfn2index(iFileId,'bzeros')
  	  iSecId = sfselect(iFileId, iSecIndex)
  	  iStatus = sfginfo(iSecId,caSec,iRank,iaDim1d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading bzero info'
  	    STOP
          ENDIF
  	  DO i=1,iRank
  	    iaStart1d(i) = 0
  	    iaEdge1d(i) = iaDim1d(i)
  	  ENDDO
  	  iStatus = sfrdata(iSecId,iaStart1d,1,iaEdge1d,azeros)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading bzero'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing bzero'
  	    STOP
          ENDIF
  
  	  ! Read in Data for timestep
  	  CALL rststr(caSec)
  	  caSec(1:10) = 'time_step_'
  	  CALL int2str(caTmp,iGetStep,iStart,iStop)
  	  caSec(11:11+iStop-iStart)=caTmp(iStart:iStop)	
  	  iSecIndex = sfn2index(iFileId,caSec)
  	  iSecId = sfselect(iFileId, iSecIndex)
  	  iAttrId = sffattr(iSecId,'time_step')
  	  iStatus = sfrnatt(iSecId,iAttrId,lstep)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading time step'
  	    STOP
          ENDIF
  	   
  	  iAttrId = sffattr(iSecId,'time')
*  	  iStatus = afgainfo(iSecId,iAttrId,'time',data_type,n_values)
  	  iStatus = sfrnatt(iSecId,iAttrId,time)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading time '
  	    STOP
          ENDIF
      stime = time
  	  WRITE(6,*) ' HDFTAKE: ',lstep,stime

  	  iStatus = sfginfo(iSecId,caSec,iRank,iaDim4d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTKAE: Failure geting step info'
  	    STOP
          ENDIF
  	  DO i=1,iRank
  	    iaStride4d(i) = 1
  	    iaStart4d(i) = 0
  	    iaEdge4d(i) = iaDim4d(i)
  	  ENDDO
  	  iStatus = sfrdata(iSecId,iaStart4d,iaStride4d,iaEdge4d,varout)
          IF (iStatus .EQ. -1) THEN
  	     WRITE(6,*) 'Failure reading step'
  	     STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing step'
  	    STOP
          ENDIF
  
  	  !Read in magnetic field data
  	  CALL rststr(caSec)
          caSec(1:13) = 'bi_time_step_'
          caSec(14:14+iStop-iStart)=caTmp(iStart:iStop)
          iSecIndex = sfn2index(iFileId,caSec)
  	  iSecId = sfselect(iFileId, iSecIndex)
          iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading bi info'
  	    STOP
          ENDIF
          DO i=1,iRank
  	    iaStride3d(i) = 1
            iaStart3d(i) = 0
            iaEdge3d(i) = iaDim3d(i)
          ENDDO
          iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,bi)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading bi'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing bi'
  	    STOP
          ENDIF
  	  CALL rststr(caSec)
          caSec(1:13) = 'bj_time_step_'
          caSec(14:14+iStop-iStart)=caTmp(iStart:iStop)
          iSecIndex = sfn2index(iFileId,caSec)
  	  iSecId = sfselect(iFileId, iSecIndex)
          iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading bj info'
  	    STOP
          ENDIF
          DO i=1,iRank
            iaStart3d(i) = 0
            iaEdge3d(i) = iaDim3d(i)
          ENDDO
          iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,bj)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading bj'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing bj'
  	    STOP
          ENDIF
  	  CALL rststr(caSec)
          caSec(1:13) = 'bk_time_step_'
          caSec(14:14+iStop-iStart)=caTmp(iStart:iStop)
          iSecIndex = sfn2index(iFileId,caSec)
  	  iSecId = sfselect(iFileId, iSecIndex)
          iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading bk info'
  	    STOP
          ENDIF
          DO i=1,iRank
            iaStart3d(i) = 0
            iaEdge3d(i) = iaDim3d(i)
          ENDDO
          iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,bk)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading bk'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing bk'
  	    STOP
          ENDIF
  
  	  !Read in electric field data
  	  CALL rststr(caSec)
          caSec(1:13) = 'ei_time_step_'
          caSec(14:14+iStop-iStart)=caTmp(iStart:iStop)
          iSecIndex = sfn2index(iFileId,caSec)
  	  iSecId = sfselect(iFileId, iSecIndex)
          iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading ei info'
  	    STOP
          ENDIF
          DO i=1,iRank
            iaStart3d(i) = 0
            iaEdge3d(i) = iaDim3d(i)
          ENDDO
          iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,etoti)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading ei'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing ei'
  	    STOP
          ENDIF
  	  CALL rststr(caSec)
          caSec(1:13) = 'ej_time_step_'
          caSec(14:14+iStop-iStart)=caTmp(iStart:iStop)
          iSecIndex = sfn2index(iFileId,caSec)
  	  iSecId = sfselect(iFileId, iSecIndex)
          iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading ej info'
  	    STOP
          ENDIF
          DO i=1,iRank
            iaStart3d(i) = 0
            iaEdge3d(i) = iaDim3d(i)
          ENDDO
          iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,etotj)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading ej'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing ej'
  	    STOP
          ENDIF
  	  CALL rststr(caSec)
          caSec(1:13) = 'ek_time_step_'
          caSec(14:14+iStop-iStart)=caTmp(iStart:iStop)
          iSecIndex = sfn2index(iFileId,caSec)
  	  iSecId = sfselect(iFileId, iSecIndex)
          iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading ek info'
  	    STOP
          ENDIF
          DO i=1,iRank
            iaStart3d(i) = 0
            iaEdge3d(i) = iaDim3d(i)
          ENDDO
          iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,etotk)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading ek'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) then
  	    WRITE(6,*) ' HDFTAKE: Failure closing ek'
  	    STOP
          ENDIF
  	
  	  ! read in ionospheric common block
   	  CALL rststr(caSec)
     	  caSec(1:14) = 'ion_time_step_'
          caSec(15:15+iStop-iStart)=caTmp(iStart:iStop)
          iSecIndex = sfn2index(iFileId,caSec)
          iSecId = sfselect(iFileId, iSecIndex)
          iStatus = sfginfo(iSecId,caSec,iRank,iaDim3d,iNtype,iNAttr)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading ionospheric info'
  	    STOP
          ENDIF
          DO i=1,iRank
            iaStart3d(i) = 0
            iaEdge3d(i) = iaDim3d(i)
          ENDDO
          iStatus = sfrdata(iSecId,iaStart3d,iaStride3d,iaEdge3d,
     $                    rionvar)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure reading ionosphere'
  	    STOP
          ENDIF
  	  iStatus = sfendacc(iSecId)
          IF (iStatus .EQ. -1) THEN
  	    WRITE(6,*) ' HDFTAKE: Failure closing ionosphere'
  	    STOP
          ENDIF
  
          DO j=1,njp1
            DO k=1,nk2p3
  	      psinorth(j,k) = psiN(j,k)
  	      curnorth(j,k) = curN(j,k)
  	      sig_p_N(j,k) =  pederN(j,k)
  	      sig_h_N(j,k) = hallN(j,k)
  	      psisouth(j,k) = psiS(j,k)
  	      cursouth(j,k) = curS(j,k)
  	      sig_p_S(j,k) = pederS(j,k)
  	      sig_h_S(j,k) = hallS(j,k)
            ENDDO
          ENDDO
    
        ENDIF

	IF (InFile) THEN
	   HDFTAKE = 0
	ELSE
	   HDFTAKE = 0
        ENDIF

	iStatus = sfend(iFileId)
        if (iStatus .EQ. -1) then
	  WRITE(6,*) ' HDFTAKE: Failure closing file'
	  STOP
        ENDIF
*	WRITE(6,*) ' HDFTAKE: Done.'

c
c Now we do some processing to create cartiesan data in the correct units
c

        CALL metric
  

c       bkress
*        in metric
*        X(I,J,K) = X2(I,J)
*        Y(I,J,K) = Y2(I,J)*cosPHI(K)
*  600   Z(I,J,K) = Y2(I,J)*sinPHI(K)
*           print *
*           do i=1,nip1
*             print *,x(i,6,8)/6378.137e5
*             print *,y(1,1,1)
*           enddo
*           stop

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

      print *
      print *,'i,j,k=',ni,1,1
      print *,'writing field values at x,y,z=',
     &X(ni,1,1)/rearth,Y(ni,1,1)/rearth,Z(ni,1,1)/rearth
      print *,'N(#/cc) is newvar(1,53,1,1)=',newvar(1,53,1,1)
      print *,'vx(km/s) is newvar(2,53,1,1)=',newvar(2,53,1,1)
      print *,'bz(nT) is newvar(8,53,1,1)=',newvar(8,53,1,1)
      print *

c
c
        if (interp) then
           print *,'calling cart1d'
           call cart1d(newvar,points,14,17,npts,allcart,intfirst)
        endif

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

c
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
