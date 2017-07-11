
      SUBROUTINE rbelt_lshell
     &(x0,t0,pa,lmonly,Lm,Lstar,leI0,bloc,bmir,bmin,err,errcode)

c rbelt L* calculator
c inputs: x,t0,pa,lmonly are position, time, pitch angle, and  
c L McIlwain only logical switch (L McIlwain calc. involves tracing the initial 
c guiding field line only. The L* calc. involves tracing the entire drift shell.)
c outputs: Lm,Lstar,leI0,bloc,bmir,bmin

      implicit none
      include 'rbelt-bounds.inc'

      logical lmonly
      integer err,errcode
      real x0(3),x(3),xmir(3),t0,rsav,pa,lm,lstar,leI0,bloc,bmir,bmin

      errcode=0

c     make sure we do not modify x0
      x(1)=x0(1)
      x(2)=x0(2)
      x(3)=x0(3)

c     temp. rmin below should be a parameter
      rsav=rmin
      rmin = 0.1

c     find mirror point
      call find_bmir(x,t0,pa,bloc,bmir,xmir,errcode)

c     get L-shell
c     original ONERA-DESP
*       SUBROUTINE calcul_Lstar(t_resol,r_resol,
*     &         lati,longi,alti,Lm,Lstar,leI0,B0,Bmin)

      if (errcode.eq.0) call calc_lshell
     &(xmir,t0,Lm,lstar,leI0,bmir,bmin,lmonly,err,errcode)
      rmin = rsav

*      if (errcode.ne.0) then
*         print *,'*** WARNING: L-star calculator failure'
*         print *,'probably an incomplete drift-shell'
*      endif

      return
      end

c MODIFIED VERSION OF ONERA-DESP L-STAR CALCULATOR -- B.KRESS, 2011
c contents of calcul_Lstar_o.f
!***************************************************************************************************
! Copyright 2003, 2004, D. Boscher, S. Bourdarie
!
! This file is part of IRBEM-LIB.
!
!    IRBEM-LIB is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    IRBEM-LIB is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
!
C Boscher modifie pour la nieme fois le 4Feb2004
C
C

       SUBROUTINE calc_lshell
     &(xsm,t0,Lm,Lstar,leI0,Bmir,Bmin,lmonly,err,errcode)

c t_resol & r_resol parameters defined below
c     t_resol -- (= options(3)+1)
c options(3rd element): resolution to compute L* (0 to 9) where 0 is the
c recomended value to ensure a good ratio precision/computation time (i.e. an
c error of ~2% at L=6). The higher the value the better will be the precision,
c the longer will be the computing time. Generally there is not much improvement
c for values larger than 4. Note that this parameter defines the integration
c step dtheta
c
c     r_resol -- (= options(4)+1) 
c coptions(4th element): resolution to compute L* (0 to 9). The higher the value
c the better will be the precision, the longer will be the computing time. It is
c recommended to use 0 (usually sufficient) unless L* is not computed on a LEO
c orbit. For LEO orbit higher values are recommended. Note that this parameter
c defines the integration step (\u03c6) along the drift shell such as
c d\u03c6=(2\u03c0)/(25*[options(4th element)+1])
c
c Takes lat, long, and altitude in geographic coordinates
c
c     lati -- second coordinate according to sysaxes. If sysaxes is 0 then
c latitude has to be in degrees otherwise use dimensionless variables (in Re)
c
c     longi -- third coordinate according to sysaxes. If sysaxes is 0 then
c longitude has to be in degrees otherwise use dimensionless variables (in Re).
c
c     alt -- ifirst coordinate according to sysaxes. If sysaxes is 0 then
c altitude has to be in km otherwise use dimensionless variables (in Re)
c
c  sysaxes: long integer to define which coordinate system is provided in
c 0: GDZ (alti, lati, East longi - km,deg,deg)
c 1: GEO (cartesian) - Re
c 2: GSM (cartesian) - Re
c 3: GSE (cartesian) - Re
c 4: SM (cartesian) - Re
c 5: GEI (cartesian) - Re
c 6: MAG (cartesian) - Re
c 7: SPH (geo in spherical) - (radial distance, lati, East longi - Re, deg., deg.)
c 8: RLL  (radial distance, lati, East longi - Re, deg., deg. - prefered than 7)
c
c     Lm
c     Lstar
c     leI0
c     B0
c     Bmin

C
       IMPLICIT NONE
       REAL     baddata
       PARAMETER  (baddata=-1.e8)
       REAL     rmin
       PARAMETER  (rmin=2.6)
C
       integer  Nreb_def,Nder_def,Ntet_def
       PARAMETER (Nreb_def = 50, Nder_def = 25, Ntet_def = 720)
C
       integer  Nder,Nreb,Ntet
       integer  n_resol,t_resol,r_resol
C
c      set resolution here
*************************************
       parameter(t_resol=1,r_resol=1)
*************************************
C
       integer  Nrebmax
       REAL     rr,rr2
       REAL     x1(3),x2(3)
       REAL     xmin(3)
C
c     bkress
       REAL     xsm(3),Bprev
C
       REAL     B(3),Bl,B0,Bmir,Bmin,B1,B3
       REAL     dsreb,smin
C
c      Ifail - get_fields failed
       integer  I,J,Iflag,Iflag_I,Ifail
       REAL     Lm,Lstar,Lb
       REAL     leI,leI0,leI1
       REAL     XY,YY
       REAL     aa,bb
C
       REAL     pi,rad,tt
       REAL     tet(10*Nder_def),phi(10*Nder_def)
       REAL     tetl,tet1,dtet
       REAL     somme
C
c      set B_zero for LM calc.
       REAL     Bo
       parameter (bo=30000.0)
C
c      bkress
       REAL t0,stetl,stetl2
       integer err,errcode
       LOGICAL lmonly
C
       Nder=Nder_def*r_resol
       Nreb=Nreb_def
       Ntet=Ntet_def*t_resol
       pi = 4.0*ATAN(1.0)
       rad = pi/180.0
       dtet = pi/Ntet
C
       Nrebmax = 20*Nreb
C
c     bkress
      err=0.
      errcode=0.
      Lm=0.
      lstar=0.
      leI0=0.
      bmir=0.
      bmin=0.
C
       rr = SQRT(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))
       tt = ACOS(xsm(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
C
C CHAMP gets total (external and internal) magnetic field
C inputs: xGEO(3)
C outputs: BxGEO(3), B0 (B mabnitude), ifail if x,y,z out of range
C

*       print *,'xsm r=',
*     & sqrt(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))

       CALL CHAMP(xsm,B,B0,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=Ifail
          err=2
	  RETURN
       ENDIF
       Bmin = B0
C
       dsreb = Lb/Nreb
C
C calcul du sens du depart
C calculation of the direction of the departure
C
*       print *,'*** 2: calculation of the direction of the departure'
       CALL sksyst (-dsreb,xsm,x1,Bl,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=Ifail
          err=2
	  RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xsm,x2,Bl,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=Ifail
          err=2
	  RETURN
       ENDIF
       B3 = Bl
C
C attention cas equatorial
C equatorial case attention
C
*       print *,'*** 3: equatorial case attention'
       IF(B1.GT.B0 .AND. B3.GT.B0)THEN
         aa = 0.5*(B3+B1-2.0*B0)
         bb = 0.5*(B3-B1)
         if (abs(-0.5*bb).gt.abs(aa)) then
            print *,'lstar: abs(-0.5*bb).gt.abs(aa)'
            stop
         endif
         smin = -0.5*bb/aa
         Bmin = B0 - aa*smin*smin
         leI0 = SQRT(1.0-Bmin/B0)*2.0*ABS(smin*dsreb)
         Lm = (Bo/Bmin)**(1.0/3.0)
         GOTO 100
       ENDIF
c      calculation of the direction of the departure here
c      go in direction of decreasing field magnitude    
       IF (B3.GT.B1) THEN
         dsreb = -dsreb
       ENDIF
C
C calcul de la ligne de champ et de I
C calculation of the line of field and I
C
*       print *,'*** 4: calculation of the line of field and I'
c      initialize bmin & leI
       Bmin = B0
       leI = 0.0
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO
C
       bprev=B0

*      print *
*      print *,'x1 r, theta, b=',
*     &sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3)),
*     &atan2(sqrt(x1(1)*x1(1)+x1(2)*x1(2)),x1(3))/rad,bl

       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)

*         print *,'x2 r, theta, b=',
*     &   sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)),
*     &   atan2(sqrt(x2(1)*x2(1)+x2(2)*x2(2)),x2(3))/rad,bl

         IF (Ifail.ne.0) THEN
            errcode=Ifail
            err=2
	    RETURN
	 ENDIF
c        retain minimum field value & position 
         IF (Bl.LT.Bmin) THEN
           xmin(1) = x2(1)
           xmin(2) = x2(2)
           xmin(3) = x2(3)
           Bmin = Bl
         ENDIF
c        if field exceeds field at mirror point then we're done 
         IF (Bl.GT.B0) GOTO 20
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
         leI = leI + SQRT(1.0-Bl/B0)
	 B1 = Bl

         Bprev=bl

       ENDDO
20     CONTINUE

C
       IF (J.GE.Nrebmax) THEN !open field line
          errcode=7
          err=2
          print *,'J.GE.Nrebmax: open field line'
          stop
	  RETURN
       ENDIF
C

*       print *,'Bprev,B1,B0,bl',Bprev,B1,B0,bl
*       print *,'SQRT(1.0-B1/B0)=',SQRT(1.0-B1/B0)
*       print *,'0.5*SQRT(1.0-B1/B0)*(B0-Bl)/(Bl-B1)=',
*     &0.5*SQRT(1.0-B1/B0)*(B0-Bl)/(Bl-B1)
*       print *,'leI,ABS(dsreb),leI0',leI,ABS(dsreb),leI*ABS(dsreb)
*       stop

c      this does not look right to me:
c      shouldn't we be adding (not subtracting) a small quantity
       leI = leI+0.5*SQRT(1.0-B1/B0)*(B0-Bl)/(Bl-B1)
       leI = leI*ABS(dsreb)
       leI0 = leI

C
C calcul de L Mc Ilwain (Mc Ilwain-Hilton)
C 
*       print *,'*** 5: calcul de L Mc Ilwain B0,B0/Bo=',B0,B0/Bo
       XY = leI*leI*leI*B0/Bo
       YY = 1.0 + 1.35047D0*XY**(1.0/3.0)
     &      + 0.465376D0*XY**(2.0/3.0)
     &      + 0.0475455D0*XY
       Lm = (Bo*YY/B0)**(1.0/3.0)
*      print *,'leI0,L McIlwain=',leI0,Lm

C
C calcul de Bmin
C
*       print *,'*** 6: calcul de Bmin'
       CALL sksyst(dsreb,xmin,x1,B3,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=Ifail
          err=2
	  RETURN
       ENDIF
       CALL sksyst(-dsreb,xmin,x1,B1,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=Ifail
          err=2
	  RETURN
       ENDIF
       aa = 0.5*(B3+B1-2.0*Bmin)
       bb = 0.5*(B3-B1)
       smin = -0.5*bb/aa
       Bmin = Bmin - aa*smin*smin
*       IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.0) THEN
c       bkress -- modified below
*        Lm = -Lm
*        print *,'conjugate point inside 1RE surface'
*       ENDIF

100    CONTINUE

c     bkress
      if (lmonly.eqv..true.) then
        IF (ABS(Lm) .GT. 10.0) then
           print *,'ABS(Lm) .GT. 10.0'
           errcode=7
           err=2
           stop
        ENDIF
        return
      ENDIF

C
C derive
C
C calcul du point sur la ligne de champ a la surface de la terre du
C cote nord
C
C calculation of the point on the line of field on the surface of the ground of
C the northern dimension
C
*       print *,'*** 7: field line point on the N. surface of the ground'
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO
       dsreb = ABS(dsreb)
       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.ne.0) THEN
            errcode=Ifail
            err=1
	    RETURN
	 ENDIF
	 rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
	 IF (rr.LT.rmin) GOTO 102
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
102    CONTINUE
       smin = sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
       smin = (rmin-smin)/(rr-smin)
       CALL sksyst(smin*dsreb,x1,x2,Bl,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=Ifail
          err=1
	  RETURN
       ENDIF
       rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
       tet(1) = ACOS(x2(3)/rr)
       phi(1) = ATAN2(x2(2),x2(1))
C
C et on tourne -> on se decale sur la surface en phi et on cherche teta
C pour avoir leI0 et B0 constants
C
C and one turns - > one shifts on surface in phi and one seeks teta to have
C leI0 and B0 constant
C
*       print *,'*** 8: find thetas around N. cap with I and Bm constan'
       dsreb = -dsreb

*       print *,'1,phi(1) =',1,phi(1)
*       print *,'1,tet(1)=',1,tet(1)/rad

       DO I = 2,Nder
        phi(I) = phi(I-1)+2.0*pi/Nder

*        print *,'i,phi(i) =',i,phi(i)

        Iflag_I = 0

c       bkress -- I do not understand original code below (commented out)
c       assumes Ilflag somewhere initialized to zero??
c       initialize/set tetl
*	IF (Ilflag.EQ.0) THEN
         tetl = tet(I-1)
	 IF (I.GT.2) tetl = 2.0*tet(I-1)-tet(I-2)
         tet1 = tetl
*	ELSE
*	 tetl = tet(I)
*	 tet1 = tetl
*	ENDIF

	leI1 = baddata
C
107     CONTINUE
c       initial point on earth's surface
        x1(1) = rmin*SIN(tetl)*COS(phi(I))
        x1(2) = rmin*SIN(tetl)*SIN(phi(I))
        x1(3) = rmin*COS(tetl)
        Iflag = 0
        leI = baddata
C
c ********* bkress ***
	    CALL CHAMP(x1,B,Bprev,Ifail,t0)
	    IF (Ifail.ne.0) THEN
               errcode=Ifail
               err=1
	       RETURN
	    ENDIF
*        print *
*        print *,'i,tetl,phi(I),Bprev=',i,tetl/rad,phi(i)/rad,Bprev
**********************
C
c       trace field line & get I
        DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.ne.0) THEN
            errcode=Ifail
            err=1
	    RETURN
	 ENDIF
         rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)

*         print *,'x2,Bl,rr',x2,Bl,sqrt(rr2)

c        when we pass the northern mirror point,
c        start calculating/integrating I
	 IF (Bl.LT.B0) THEN
c         initialize I (first time only, Iflag initially 0, set to 1 below)
	  IF (Iflag .EQ. 0) THEN
	    CALL CHAMP(x1,B,B1,Ifail,t0)
	    IF (Ifail.ne.0) THEN
               errcode=Ifail
               err=1
	       RETURN
	    ENDIF
	    leI = 0.5*SQRT(1.0-Bl/B0)*(1.0+(Bl-B0)/(Bl-B1))
	    Iflag = 1
	  ELSE
	    leI = leI+SQRT(1.0-Bl/B0)
	  ENDIF
	 ENDIF
c        break after we pass southern mirror point
         IF (Bl.GT.B0 .AND. Iflag.EQ.1) GOTO 103
	 IF (rr2.LT.rmin) GOTO 103
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
        ENDDO
103     CONTINUE

*        print *,'B0,Bl,leI0,leI =',B0,Bl,leI0,leI

c Pourquoi?
c       stepped inside of r=1 to find B > B0 (to get to southern mirror point)
c       or B reached southern cap and B never went below B0 (theta too large --
c       but close enough? we use this theta.)
        IF (rr2.LT.rmin) THEN
           leI = baddata
	ENDIF

c       bkress -- modified the following
        IF (J.LT.Nrebmax .AND. rr2.GE.rmin) THEN
            CALL CHAMP(x1,B,B1,Ifail,t0)
	    IF (Ifail.ne.0) THEN
               errcode=Ifail
               err=1
	       RETURN
	    ENDIF
            leI = leI+0.5*SQRT(1.0-B1/B0)*(B0-Bl)/(Bl-B1)
            leI = leI*ABS(dsreb)
	ENDIF

c       for first tetl only
c       (Iflag_I set to 1 below)
        IF (Iflag_I .EQ.0) THEN
c        ?? if we did not find southern mirror pt. (B0) ??
	 IF (J.GE.Nrebmax) THEN
	  tetl = tetl-dtet
         ELSE
	  tetl = tetl+dtet
	 ENDIF
	 leI1 = leI
	 tet1 = tetl
	 Iflag_I = 1
	 GOTO 107
	ENDIF

c       if leI0 is between leI and leI1,
c       we are done looking for correct theta
c       (note, leI1 set = leI above,
c       we do not goto 108 on the first time through)
	IF ((leI-leI0)*(leI1-leI0) .LT. 0.0) GOTO 108

	leI1 = leI
	tet1 = tetl
c       *** set new tetl here ***
	IF (leI.LT.leI0) THEN
c        if leI too small, move poleward
	 tetl = tetl-dtet
	ElSE
c        if leI too big, move equatorward
	 tetl = tetl+dtet
	ENDIF

	IF (tetl.GT.pi .OR. tetl.LT.0.0) GOTO 108
	GOTO 107
108     CONTINUE

        tet(I) = 0.5*(tetl+tet1)

c       done getting correct theta **********************
*        print *,'i,tet(I)=',i,tet(I)/rad

	IF (J.GE.Nrebmax .AND. leI.GT.0.0) THEN
         errcode=7
         err=1
*	  print *,'J.GE.Nrebmax .AND. leI.GT.0.0'
*         stop
	 RETURN
	ENDIF
C
        x1(1) = rmin*SIN(tet(I))*COS(phi(I))
        x1(2) = rmin*SIN(tet(I))*SIN(phi(I))
        x1(3) = rmin*COS(tet(I))

	CALL CHAMP(x1,B,Bl,Ifail,t0)
	IF (Ifail.ne.0) THEN
           errcode=Ifail
           err=1
	   RETURN
	ENDIF
	
*	print *,'Bl=',Bl
	
c       because of the method choosen (i.e. look for correct theta as we 
c       go around the polar cap at 1 RE) we have to restrict to L-shells 
c       that stay above 1RE. Here we check to see if the northern mirror
c       point is inside of 1RE.
        IF (Bl.LT.B0) THEN
         errcode=7
         err=1
*	 print *,'Bl.LT.B0',Bl,B0
*	 print *,'xsm=',xsm
*        stop
	 RETURN
        ENDIF
       ENDDO
C
C calcul de somme de BdS sur la calotte nord
C calculation of nap of BdS on the northern cap
C
*       print *,'*** 9: calc. sum of BdS on the northern cap'

       x1(1) = 0.0
       x1(2) = 0.0
       x1(3) = rmin
       CALL CHAMP(x1,B,Bl,Ifail,t0)
       IF (Ifail.ne.0)THEN
           errcode=Ifail
           err=1
	   RETURN
       ENDIF
       somme = abs(x1(1)*B(1)+x1(2)*B(2)+x1(3)*B(3))*pi*rmin*rmin*
     & dtet*dtet/4.0
       DO I = 1,Nder
         tetl = 0.0
         DO J = 1,Ntet
	  tetl = tetl+dtet
	  IF (tetl .GT. tet(I)) GOTO 111
          x1(1) = rmin*SIN(tetl)*COS(phi(I))
          x1(2) = rmin*SIN(tetl)*SIN(phi(I))
          x1(3) = rmin*COS(tetl)
          CALL CHAMP(x1,B,Bl,Ifail,t0)
          IF (Ifail.ne.0)THEN
              errcode=Ifail
              err=1
	      RETURN
          ENDIF
          somme = somme + abs((x1(1)*B(1)+x1(2)*B(2)+x1(3)*B(3)))*
     &    rmin*SIN(tetl)*dtet*2.D0*pi/Nder
	 ENDDO
111      CONTINUE
       ENDDO
       Lstar = 2.0*pi*Bo/somme

       return
       END

C***********************************************************************

      SUBROUTINE CHAMP(y,bxsm,Bl,Ifail,t)
C
      IMPLICIT NONE
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-status.inc'
      include 'rbelt-grid.inc'
      integer   Ifail
      REAL      y(3),Bxsm(3),Bl,t

*      print *
*      print *,'in SUBROUTINE CHAMP'

c     initialize ifail
      Ifail=0

      if (sys.eq.2) then
         print *,'need to rotate to geo coordinates'
         stop
      endif

      status=0
      call get_fields(y,t)
      if (status.eq.1 .or. status.eq.2) then
*         print *,'hit inner or outer boundary in L-shell calc.',status
*         print *,'x,y,z,r=',
*     &   y(1),y(2),y(3),sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
*         stop
         Ifail=status
         return
      endif

c     b in sm coordinate system
      Bxsm(1)=bx/ffactor/ntg
      Bxsm(2)=by/ffactor/ntg
      Bxsm(3)=bz/ffactor/ntg
      Bl=sqrt(bxsm(1)*bXsm(1)+bxsm(2)*bxsm(2)+bxsm(3)*bxsm(3))

c      if (sys.eq.2) print *,'need to rotate B to geo coordinates'

      return
      end

C***********************************************************************

      SUBROUTINE CHAMP2(y,bxsm,Bl,Ifail,t)
C
      IMPLICIT NONE
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-status.inc'
      include 'rbelt-grid.inc'
      integer   Ifail
      REAL      y(3),Bxsm(3),Bl,t

*      print *
*      print *,'in SUBROUTINE CHAMP'

c     initialize ifail
      Ifail=0

      if (sys.eq.2) then
         print *,'need to rotate to geo coordinates'
         stop
      endif

      status=0
      call get_fields2(y,t)
      if (status.eq.1 .or. status.eq.2) then
*         print *,'hit inner or outer boundary in L-shell calc.',status
*         print *,'x,y,z,r=',
*     &   y(1),y(2),y(3),sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
*         stop
         Ifail=status
         return
      endif

c     b in sm coordinate system
      Bxsm(1)=bx/ffactor/ntg
      Bxsm(2)=by/ffactor/ntg
      Bxsm(3)=bz/ffactor/ntg
      Bl=sqrt(bxsm(1)*bXsm(1)+bxsm(2)*bxsm(2)+bxsm(3)*bxsm(3))

      return
      end

C
C***********************************************************************
C* RUNGE-KUTTA d'ordre 4
C***********************************************************************
        SUBROUTINE sksyst(h,xx,x2,Bl,Ifail,t0)
C
        IMPLICIT NONE
C
        integer  Ifail
        REAL xx(3),x2(3)
        REAL B(3),Bl
        REAL h
        REAL xwrk(4,3)
        REAL t0
C
C-----------------------------------------------------------------------
C
c        write(6,*)'sksyst'
c        write(6,*)xx(1),xx(2),xx(3),h
        CALL CHAMP(xx,B,Bl,Ifail,t0)
c	write(6,*)xx,B,Bl,Ifail
	IF (Ifail.ne.0) RETURN
c        write(6,*)'b',B(1),B(2),B(3),Bl
        xwrk(1,1) = h*B(1)/Bl
        xwrk(1,2) = h*B(2)/Bl
        xwrk(1,3) = h*B(3)/Bl
        x2(1) = xx(1)+xwrk(1,1)/2.0
        x2(2) = xx(2)+xwrk(1,2)/2.0
        x2(3) = xx(3)+xwrk(1,3)/2.0
c        write(6,*)x2(1),x2(2),x2(3),Bl
C
        CALL CHAMP(x2,B,Bl,Ifail,t0)
	IF (Ifail.ne.0) RETURN
        xwrk(2,1) = h*B(1)/Bl
        xwrk(2,2) = h*B(2)/Bl
        xwrk(2,3) = h*B(3)/Bl
        x2(1) = xx(1)+xwrk(2,1)/2.0
        x2(2) = xx(2)+xwrk(2,2)/2.0
        x2(3) = xx(3)+xwrk(2,3)/2.0
c        write(6,*)x2(1),x2(2),x2(3),Bl
C
        CALL CHAMP(x2,B,Bl,Ifail,t0)
	IF (Ifail.ne.0) RETURN
        xwrk(3,1) = h*B(1)/Bl
        xwrk(3,2) = h*B(2)/Bl
        xwrk(3,3) = h*B(3)/Bl
        x2(1) = xx(1)+xwrk(3,1)
        x2(2) = xx(2)+xwrk(3,2)
        x2(3) = xx(3)+xwrk(3,3)
c        write(6,*)x2(1),x2(2),x2(3),Bl
C
        CALL CHAMP(x2,B,Bl,Ifail,t0)
	IF (Ifail.ne.0) RETURN
        xwrk(4,1) = h*B(1)/Bl
        xwrk(4,2) = h*B(2)/Bl
        xwrk(4,3) = h*B(3)/Bl
C
C-----------------------------------------------------------------------
C
        x2(1) = xx(1)+(   xwrk(1,1)+2.0*xwrk(2,1)
     &               + 2.0*xwrk(3,1)+   xwrk(4,1))/6.0
        x2(2) = xx(2)+(   xwrk(1,2)+2.0*xwrk(2,2)
     &               + 2.0*xwrk(3,2)+   xwrk(4,2))/6.0
        x2(3) = xx(3)+(   xwrk(1,3)+2.0*xwrk(2,3)
     &               + 2.0*xwrk(3,3)+   xwrk(4,3))/6.0
        CALL CHAMP(x2,B,Bl,Ifail,t0)
	IF (Ifail.ne.0) RETURN
c        write(6,*)x2(1),x2(2),x2(3),Bl
c        read(5,*)
C
        RETURN
        END

!***************************************************************************************************
! Copyright 2004 S. Bourdarie
!
! This file is part of IRBEM-LIB.
!
!    IRBEM-LIB is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    IRBEM-LIB is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
!
C S. Bourdarie (June 2004)
c Modified S. Bourdarie (July 2004)
c Modified bkress (2011-2013)

       SUBROUTINE find_bmir(xsm,t0,alpha,Bposit,Bmir,xmir,errcode)
C
       IMPLICIT NONE
C
       integer  Nreb
       PARAMETER (Nreb = 50)
C
       integer  Ifail,errcode
       integer  Nrebmax
       real     rr
       real     xsm(3),t0,x1(3),x2(3)
       real     xmir(3)
*       real     B1,B3
       real     B(3),Bl,B0
       real     dsreb
       integer  I,J
       real     Lb
       real     pi,rad
       real     tt
       real     cste
       real     Bposit,Bmir
       real     sn2,alpha
C
       pi = 4.0*ATAN(1.0)
       rad = pi/180.0
C
       Nrebmax = 20*Nreb
C
       rr = SQRT(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))
       tt = ACOS(xsm(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
C
       errcode=0
C
       CALL CHAMP(xsm,B,B0,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=ifail
	  RETURN
       ENDIF
       Bmir=B0
       Bposit=B0
       sn2=SIN(alpha)*SIN(alpha)
       cste=B0/sn2

c      set field line tracing step size
       dsreb = Lb/(Nreb*1.0)

c      if alpha=90 deg. we are done
       if (alpha.eq.pi/2) then
         do I = 1,3
           xmir(I)  = xsm(I)
         enddo
         return
c      trace in direction specified by PA
       elseif (alpha.gt.pi/2.) then
         dsreb = -dsreb
       endif

C calcul de la ligne de champ
C (calculation of the line of field)

c      trace field line in direction of parallel velocity
c      and set x1 = position preceding B/B_mirror > 1
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO
       DO J = 1,Nrebmax
c        x1 unchanged by sksyst, x2 gets new/advanced position
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.ne.0) THEN
            errcode=ifail
	    RETURN
         ENDIF
c        this is current Bl/B_mirror
	 sn2=Bl/cste
c        exit loop when B/B_mirror > 1
c        x1 is position immediately preceding
*         print *,'   x2,Bl,sn2=',x2,Bl,sn2
	 if (sn2.GT.1.0) GOTO 20
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
20     CONTINUE
C
       IF (J.GE.Nrebmax) THEN !open field line
          errcode=2
	  RETURN
       ENDIF

C           
C calcul de Bmirror point
C

c      now approach mirror point with smaller steps
       DO I = 1,3
         xmir(I)  = x1(I)
       ENDDO
       CALL CHAMP(x1,B,Bl,Ifail,t0)
       DO i=1,100
         CALL sksyst(dsreb/100.0,xmir,x2,Bl,Ifail,t0)
         IF (Ifail.ne.0) THEN
            errcode=ifail
	    RETURN
         ENDIF
	 sn2=Bl/cste
c        again stop with xmir immediately preceding B/B_mirror > 1 location
	 if (sn2 .GT.1.0) GOTO 30
         xmir(1) = x2(1)
         xmir(2) = x2(2)
         xmir(3) = x2(3)
	 Bmir=Bl
       ENDDO
30     CONTINUE
c      that's all we need

       RETURN
       END

C
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
       SUBROUTINE find_bmin(xsm,t0,Bloc,Bmin,xmin,errcode)
C
       IMPLICIT NONE
C
       integer  Nreb
       PARAMETER (Nreb = 50)
C
       integer  Ifail,errcode
       integer  Nrebmax
       real     rr
       real     xsm(3),t0,x0(3),x1(3),x2(3)
       real     xmin(3)
       real     B(3),Bl,B0,B1,B3
       real     dsreb
       integer  I,J,ind
       real     Lb
       real     pi,rad
       real     tt
       real     Bloc,Bmin
       REAL     aa,bb,smin
C
       pi = 4.0*ATAN(1.0)
       rad = pi/180.0
C
       Nrebmax = 20*Nreb
C
       rr = SQRT(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))
       tt = ACOS(xsm(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
C
       errcode=0
C
       CALL CHAMP(xsm,B,B0,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=ifail
	  RETURN
       ENDIF
       Bloc=B0

c      set field line tracing step size
       dsreb = Lb/(Nreb*1.0)

C
C calcul du sens du depart 
C (calculation of the direction of the departure)
C

       CALL sksyst(-dsreb,xsm,x1,Bl,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=ifail
	  RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xsm,x2,Bl,Ifail,t0)
       IF (Ifail.ne.0) THEN
          errcode=ifail
	  RETURN
       ENDIF
       B3 = Bl

C
C attention cas equatorial
C equatorial case attention
C
       IF(B1.GT.B0 .AND. B3.GT.B0)THEN
         aa = 0.5*(B3+B1-2.0*B0)
         bb = 0.5*(B3-B1)
         if (abs(-0.5*bb).gt.abs(aa)) then
            print *,'lstar: abs(-0.5*bb).gt.abs(aa)'
            stop
         endif
         smin = -0.5*bb/aa
         Bmin = B0 - aa*smin*smin
         xmin(1) = xsm(1) + smin*(x2(1)-x1(1))/2.
         xmin(2) = xsm(2) + smin*(x2(2)-x1(2))/2.
         xmin(3) = xsm(3) + smin*(x2(3)-x1(3))/2.
         return
       ENDIF

c      calculation of the direction of the departure here
c      go in direction of decreasing field magnitude    
       IF (B3.GT.B1) THEN
         dsreb = -dsreb
       ENDIF

C
C calcul de la ligne de champ et de I
C calculation of the line of field and I
C

c      initialize bmin
       Bmin = B0
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO

       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.ne.0) THEN
            errcode=ifail
	    RETURN
	 ENDIF
c        retain position immediately preceding bmin location
c        also, bin & bmin location
         IF (Bl.LT.Bmin) THEN
           x0(1) = x1(1)
           x0(2) = x1(2)
           x0(3) = x1(3)
           xmin(1) = x2(1)
           xmin(2) = x2(2)
           xmin(3) = x2(3)
           Bmin = Bl
         ENDIF
c        if field exceeds field at mirror point then we're done
c        could include test for double minimum here
         IF (Bl.GT.B0) GOTO 20
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
20     CONTINUE
       IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.0) THEN
        print *,'conjugate point inside 1RE surface'
        errcode=1
       ENDIF
       IF (J.GE.Nrebmax) THEN !open field line
          errcode=2
	  RETURN
       ENDIF

c      now approach & hopfully pass bmin with smaller steps starting at x0
       DO i=1,200
         CALL sksyst(dsreb/100.0,x0,x1,Bl,Ifail,t0)
         IF (Ifail.ne.0) THEN
            errcode=ifail
	    RETURN
         ENDIF
c        retain minimum field value & position 
         IF (Bl.LT.Bmin) THEN
           xmin(1) = x1(1)
           xmin(2) = x1(2)
           xmin(3) = x1(3)
           Bmin = Bl
         ENDIF
	 x0(1) = x1(1)
	 x0(2) = x1(2)
	 x0(3) = x1(3)
       ENDDO

       return
       END

C
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
       SUBROUTINE gbdbzero(xsm,t0,gbdbloc,gbdb0,x0,ifail)
C
       IMPLICIT NONE
       include 'rbelt-fields.inc'
       real gbdb,gbdb0,gbdb1,gbdb2,gbdb_prev
C
       integer  Nreb
       PARAMETER (Nreb = 50)
C
       integer  Ifail
       integer  Nrebmax
       integer  I,J,ind
       real     rr
       real     xsm(3),t0,x0(3),x1(3),x2(3)
       real     B3(3),Bl,B1,B2
       real     dsreb
       real     Lb
       real     pi,rad
       real     tt
       real     gbdbloc
       real     delta
C
       pi = 4.0*ATAN(1.0)
       rad = pi/180.0
C
       Nrebmax = 20*Nreb
C
       rr = SQRT(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))
       tt = ACOS(xsm(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
C
       CALL CHAMP2(xsm,B3,Bl,Ifail,t0)
       IF (Ifail.ne.0) RETURN
       gbdbloc=dbdx*bx+dbdy*by+dbdz*bz

c      set field line tracing step size
       dsreb = Lb/(Nreb*1.0)

C
C calcul du sens du depart 
C (calculation of the direction of the departure)
C

       CALL sksyst(dsreb,xsm,x2,B2,Ifail,t0)
       CALL CHAMP2(x2,B3,B2,Ifail,t0)
       IF (Ifail.ne.0) RETURN
       gbdb2=dbdx*bx+dbdy*by+dbdz*bz

*       print *
*       print *,'xsm',xsm
*       print *,'gbdbloc,gbdb2=',gbdbloc,gbdb2
*       print *,'Bl,B2=',Bl,B2
*       print *

       IF ((abs(gbdbloc).LT.abs(gbdb2)).and.(gbdbloc*gbdb2.gt.0.)) THEN
         dsreb = -dsreb
*         print *,'reverse direction'
       ENDIF

*       stop

C
C calcul de la ligne de champ et de I
C calculation of the line of field and I
C

c      initialize gbdb
       gbdb_prev = gbdbloc
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO

       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         CALL CHAMP2(x2,B3,Bl,Ifail,t0)
         IF (Ifail.ne.0) RETURN
         gbdb=dbdx*bx+dbdy*by+dbdz*bz

*         print *,'*** x2,B,bgbdb=',x2,Bl,gbdb

c        if gbdb goes through zero
         IF ((gbdb*gbdb_prev.le.0.).and.(gbdb.ne.0.)) GOTO 20
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
         gbdb_prev = gbdb
       ENDDO
       Ifail=7
       RETURN
20     CONTINUE

       IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.0) THEN
        print *,'conjugate point inside 1RE surface'
       ENDIF
       IF (J.GE.Nrebmax) RETURN

c      now approach & hopfully pass gbdb0 with smaller steps starting at x0
*       print *
       dsreb=dsreb/100.
       DO i=1,200
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         CALL CHAMP2(x2,B3,Bl,Ifail,t0)
         IF (Ifail.ne.0) RETURN
         gbdb=dbdx*bx+dbdy*by+dbdz*bz

*         print *
*         print *,'x2,gbdb=',x2,gbdb
*         print *,'bx,by,bz,dbdx,dbdy,dbdz=',bx,by,bz,dbdx,dbdy,dbdz

c        if gbdb goes through zero
         IF ((gbdb*gbdb_prev.le.0.).and.(gbdb.ne.0.)) GOTO 30
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
         gbdb_prev = gbdb
       ENDDO
       Ifail=7
       RETURN
30     CONTINUE

       delta=gbdb_prev/(gbdb_prev-gbdb)
       gbdb0=gbdb_prev+delta*(gbdb-gbdb_prev)
       x0(1) = x1(1) + delta*(x2(1)-x1(1))
       x0(2) = x1(2) + delta*(x2(2)-x1(2))
       x0(3) = x1(3) + delta*(x2(3)-x1(3))

       return
       END

C
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
       SUBROUTINE trace2eqtr(t0,xsm,Beqtr,ifail)
C
       IMPLICIT NONE
       include 'rbelt-fields.inc'
C
       integer  Nreb
       PARAMETER (Nreb = 50)
C
       integer  Ifail
       integer  Nrebmax
       real     rr
       real     xsm(3),t0,x0(3),x1(3),x2(3)
       real     Bl,B1,Beqtr
       real     dsreb
       integer  I,J,ind
       real     Lb
       real     pi,rad
       real     tt
       REAL     delta
C
       ifail=0
C
       pi = 4.0*ATAN(1.0)
       rad = pi/180.0
C
       Nrebmax = 20*Nreb
C
C      set field line tracing step size
       rr = SQRT(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))
       tt = ACOS(xsm(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
       dsreb = Lb/(Nreb*1.0)

       IF (xsm(3).GT.0.0) THEN
         dsreb = -dsreb
       ENDIF

C
C calcul de la ligne de champ et de I
C calculation of the line of field and I
C
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO

       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.ne.0) RETURN
c        if we pass equatorial plane
         IF ((x2(3)*x1(3).le.0.).and.(x2(3).ne.0.)) GOTO 20
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
       Ifail=7
       RETURN
20     CONTINUE

       IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.0) THEN
          print *,'conjugate point inside 1RE surface'
          stop
       ENDIF
       IF (J.GE.Nrebmax) RETURN

c      now approach & hopfully pass z=0 with smaller steps starting at x0
       b1=bl
       dsreb=dsreb/100.
       DO i=1,200
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.ne.0) RETURN
c        if z goes through zero
         IF ((x2(3)*x1(3).le.0.).and.(x2(3).ne.0.)) GOTO 30
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
         b1=bl
       ENDDO
       Ifail=7
       RETURN
30     CONTINUE

       delta=x1(3)/(x1(3)-x2(3))
       xsm(1) = x1(1) + delta*(x2(1)-x1(1))
       xsm(2) = x1(2) + delta*(x2(2)-x1(2))
       xsm(3) = x1(3) + delta*(x2(3)-x1(3))
       beqtr=b1+delta*(bl-b1)

       return
       END





