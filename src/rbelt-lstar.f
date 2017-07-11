
      SUBROUTINE rbelt_lshell(x,t0,pa,Lm,Lstar,leI0,bloc,bmir,bmin)

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-flag.inc'

      logical lmonly
      integer lmfail,lsfail,ifail
      real*8 x(3),xmir(3),t0,rmin_tmp,pa,lm,lstar,leI0,bloc,bmir,bmin

      lmonly=.false.
      if (lcalc.eq.1) lmonly=.true.

      rmin_tmp=rmin
      rmin = 0.25

      call find_bmir(x,t0,pa,bloc,bmir,xmir,ifail)

c     original ONERA-DESP
*       SUBROUTINE calcul_Lstar(t_resol,r_resol,
*     &         lati,longi,alti,Lm,Lstar,leI0,B0,Bmin)
      call calc_lshell
     &(xmir,t0,Lm,lstar,leI0,bmir,bmin,lmonly,lmfail,lsfail)
      rmin = rmin_tmp

      if ((lmfail.eq.1).or.(lsfail.eq.1)) then
         print *,'*** WARNING: L-star calculator failure'
         print *,'probably incomplete drift-shell'
         print *,'lmfail,lsfail=',lmfail,lsfail
*         stop
      else
*        print *,'leI0, L-star=',leI0,lstar
      endif

      return
      end

! MODIFIED VERSION OF ONERA-DESP L-STAR CALCULATOR -- BKRESS, 2011
! contents of calcul_Lstar_o.f
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
! Boscher modifie pour la nieme fois le 4Feb2004
!
!

       SUBROUTINE calc_lshell
     &(xsm,t0,Lm,Lstar,leI0,B0,Bmin,lmonly,lmfail,lsfail)

! t_resol & r_resol parameters defined below
!     t_resol -- (= options(3)+1)
! options(3rd element): resolution to compute L* (0 to 9) where 0 is the
! recomended value to ensure a good ratio precision/computation time (i.e. an
! error of ~2% at L=6). The higher the value the better will be the precision,
! the longer will be the computing time. Generally there is not much improvement
! for values larger than 4. Note that this parameter defines the integration
! step dtheta
!
!     r_resol -- (= options(4)+1) 
! coptions(4th element): resolution to compute L* (0 to 9). The higher the value
! the better will be the precision, the longer will be the computing time. It is
! recommended to use 0 (usually sufficient) unless L* is not computed on a LEO
! orbit. For LEO orbit higher values are recommended. Note that this parameter
! defines the integration step (\u03c6) along the drift shell such as
! d\u03c6=(2\u03c0)/(25*[options(4th element)+1])

! Takes lat, long, and altitude in geographic coordinates
!
!     lati -- second coordinate according to sysaxes. If sysaxes is 0 then
! latitude has to be in degrees otherwise use dimensionless variables (in Re)
!
!     longi -- third coordinate according to sysaxes. If sysaxes is 0 then
! longitude has to be in degrees otherwise use dimensionless variables (in Re).
!
!     alt -- ifirst coordinate according to sysaxes. If sysaxes is 0 then
! altitude has to be in km otherwise use dimensionless variables (in Re)
!
!  sysaxes: long integer to define which coordinate system is provided in
!c 0: GDZ (alti, lati, East longi - km,deg,deg)
!c 1: GEO (cartesian) - Re
!c 2: GSM (cartesian) - Re
!c 3: GSE (cartesian) - Re
!c 4: SM (cartesian) - Re
!c 5: GEI (cartesian) - Re
!c 6: MAG (cartesian) - Re
!c 7: SPH (geo in spherical) - (radial distance, lati, East longi - Re, deg., deg.)
!c 8: RLL  (radial distance, lati, East longi - Re, deg., deg. - prefered than 7)
!c
!c     Lm
!c     Lstar
!c     leI0
!c     B0
!c     Bmin

!C
       IMPLICIT NONE
       REAL*8     baddata
       PARAMETER  (baddata=-1.d31)
!C
       integer  Nreb_def,Nder_def,Ntet_def
       PARAMETER (Nreb_def = 50, Nder_def = 25, Ntet_def = 720)
!C
       integer  Nder,Nreb,Ntet
       integer  n_resol,t_resol,r_resol
!C
!c      set resolution here
!*************************************
       parameter(t_resol=1,r_resol=1)
!*************************************
!C
       integer  Nrebmax
       REAL*8     rr,rr2
       REAL*8     x1(3),x2(3)
       REAL*8     xmin(3)
!C
!c     bkress
       REAL*8     xsm(3),Bprev
!C
       REAL*8     B(3),Bl,B0,Bmin,B1,B3
       REAL*8     dsreb,smin
!C
!c      Ifail - get_fields failed
!c      Ilflag - L-shell calc. failed
       integer  I,J,Iflag,Iflag_I,Ilflag,Ifail
       REAL*8     Lm,Lstar,Lb
       REAL*8     leI,leI0,leI1
       REAL*8     XY,YY
       REAL*8     aa,bb
!C
       REAL*8     pi,rad,tt
       REAL*8     tet(10*Nder_def),phi(10*Nder_def)
       REAL*8     tetl,tet1,dtet
       REAL*8     somme
!C
!c      set B_zero for LM calc.
       REAL*8     Bo
       parameter (bo=30000.0)
!C
!c      bkress
       REAL*8 t0,stetl,stetl2
       integer lmfail,lsfail
       LOGICAL lmonly
!C
       Nder=Nder_def*r_resol
       Nreb=Nreb_def
       Ntet=Ntet_def*t_resol
       pi = 4.0*ATAN(1.0)
       rad = pi/180.0
       dtet = pi/Ntet
!C
       Nrebmax = 20*Nreb
!C
       Lm = baddata
       Lstar = baddata
       leI0 = baddata
!C
!c     bkress
      lmfail=0
      lsfail=0
!C
       rr = DSQRT(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))
       tt = DACOS(xsm(3)/rr)
       Lb  = rr/DSIN(tt)/DSIN(tt)
!c       write(6,*)'L bete ',Lb
!C
!C CHAMP gets total (external and internal) magnetic field
!C inputs: xGEO(3)
!C outputs: BxGEO(3), B0 (B mabnitude), ifail if x,y,z out of range
!C
       CALL CHAMP(xsm,B,B0,Ifail,t0)
       IF (Ifail.LT.0) THEN
          B0=baddata
          leI0 = baddata
          Bmin = baddata
          Ilflag = 0
          lmfail=1
          lsfail=1
          RETURN
       ENDIF
       Bmin = B0
!C
       dsreb = Lb/Nreb
!C
!C calcul du sens du depart
!C calculation of the direction of the departure
!C
!      print *,'*** 2: calculation of the direction of the departure'
       CALL sksyst (-dsreb,xsm,x1,Bl,Ifail,t0)
       IF (Ifail.LT.0) THEN
          leI0 = baddata
          Bmin = baddata
          Ilflag = 0
          lmfail=1
          lsfail=1
          RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xsm,x2,Bl,Ifail,t0)
       IF (Ifail.LT.0) THEN
          leI0 = baddata
          Bmin = baddata
          Ilflag = 0
          lmfail=1
          lsfail=1
          RETURN
       ENDIF
       B3 = Bl
!C
!C attention cas equatorial
!C equatorial case attention
!C
!      print *,'*** 3: equatorial case attention'
       IF(B1.GT.B0 .AND. B3.GT.B0)THEN
         aa = 0.5*(B3+B1-2.0*B0)
         bb = 0.5*(B3-B1)
         if (dabs(-0.5*bb).gt.dabs(aa)) then
            print *,'lstar: abs(-0.5*bb).gt.abs(aa)'
            stop
         endif
         smin = -0.5*bb/aa
         Bmin = B0 - aa*smin*smin
         leI0 = DSQRT(1.0-Bmin/B0)*2.0*DABS(smin*dsreb)
         Lm = (Bo/Bmin)**(1.0/3.0)
!c         write(6,*)'L McIlwain eq ',B0,leI0,Lm

!c        bkress -- only need Lm for penumbra code
         if (lmonly.eqv..true.) return

         GOTO 100
       ENDIF
!c      calculation of the direction of the departure here
!c      go in direction of decreasing field magnitude    
       IF (B3.GT.B1) THEN
         dsreb = -dsreb
       ENDIF
!C
!C calcul de la ligne de champ et de I
!C calculation of the line of field and I
!C
!      print *,'*** 4: calculation of the line of field and I'
!c      initialize bmin & leI
       Bmin = B0
       leI = 0.0
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO
!C
       bprev=B0

       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.LT.0) THEN
            leI0 = baddata
            Bmin = baddata
            Ilflag = 0
            lmfail=1
            lsfail=1
            print *,'Ifail.LT.0'
            RETURN
         ENDIF
!c        write(6,*)J,x1(1),x1(2),x1(3),Bl
!c        retain minimum field value & position 
         IF (Bl.LT.Bmin) THEN
           xmin(1) = x2(1)
           xmin(2) = x2(2)
           xmin(3) = x2(3)
           Bmin = Bl
         ENDIF
!c        if field exceeds field at mirror point then we're done 
         IF (Bl.GT.B0) GOTO 20
         x1(1) = x2(1)
         x1(2) = x2(2)
         x1(3) = x2(3)
         leI = leI + DSQRT(1.0-Bl/B0)
         B1 = Bl

         Bprev=bl

       ENDDO
20     CONTINUE

!C
       IF (J.GE.Nrebmax) THEN !open field line
          leI0 = baddata
          Bmin = baddata
          Ilflag = 0
          lmfail=1
          lsfail=1
          print *,'J.GE.Nrebmax: open field line'
          RETURN
       ENDIF
!C

!*       print *,'Bprev,B1,B0,bl',Bprev,B1,B0,bl
!*       print *,'SQRT(1.0-B1/B0)=',SQRT(1.0-B1/B0)
!*       print *,'0.5*SQRT(1.0-B1/B0)*(B0-Bl)/(Bl-B1)=',
!*     &0.5*SQRT(1.0-B1/B0)*(B0-Bl)/(Bl-B1)
!*       print *,'leI,ABS(dsreb),leI0',leI,ABS(dsreb),leI*ABS(dsreb)
!*       stop

!c      this does not look right to me:
!c      shouldn't we be adding (not subtracting) a small quantity
       leI = leI+0.5*DSQRT(1.0-B1/B0)*(B0-Bl)/(Bl-B1)
       leI = leI*DABS(dsreb)
       leI0 = leI

!C
!C calcul de L Mc Ilwain (Mc Ilwain-Hilton)
!C 
!      print *,'*** 5: calcul de L Mc Ilwain B0,B0/Bo=',B0,B0/Bo
       XY = leI*leI*leI*B0/Bo
       YY = 1.0 + 1.35047D0*XY**(1.0/3.0)
     &      + 0.465376D0*XY**(2.0/3.0)
     &      + 0.0475455D0*XY
       Lm = (Bo*YY/B0)**(1.0/3.0)
!c       write(6,*)'L McIlwain ',B0,leI0,Lm
!*      print *,'leI0,L McIlwain=',leI0,Lm

!C
!C calcul de Bmin
!C
!      print *,'*** 6: calcul de Bmin'
       CALL sksyst(dsreb,xmin,x1,B3,Ifail,t0)
       IF (Ifail.LT.0) THEN
          Bmin = baddata
          Ilflag = 0
          lsfail = 1
          RETURN
       ENDIF
       CALL sksyst(-dsreb,xmin,x1,B1,Ifail,t0)
       IF (Ifail.LT.0) THEN
          Bmin = baddata
          Ilflag = 0
          lsfail = 1
          RETURN
       ENDIF
       aa = 0.5*(B3+B1-2.0*Bmin)
       bb = 0.5*(B3-B1)
       smin = -0.5*bb/aa
       Bmin = Bmin - aa*smin*smin
       IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.0) THEN
!c       bkress -- modified below
!*        Lm = -Lm
        print *,'conjugate point inside 1RE surface'
       ENDIF

!c     bkress -- only need Lm for penumbra code
      if (lmonly.eqv..true.) return

!C
100    CONTINUE

       IF (DABS(Lm) .GT. 10.0) THEN
        Ilflag = 0
        lsfail = 1
        RETURN
       ENDIF

!C
!C derive
!C
!C calcul du point sur la ligne de champ a la surface de la terre du
!C cote nord
!C
!C calculation of the point on the line of field on the surface of the ground of
!C the northern dimension
!C
!      print *,'*** 7: field line point on the N. surface of the ground'
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO
       dsreb = DABS(dsreb)
       DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
            lsfail = 1
            RETURN
         ENDIF
         rr = dsqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
         IF (rr.LT.1.0) GOTO 102
         x1(1) = x2(1)
         x1(2) = x2(2)
         x1(3) = x2(3)
       ENDDO
102    CONTINUE
       smin = dsqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
       smin = (1.0-smin)/(rr-smin)
       CALL sksyst(smin*dsreb,x1,x2,Bl,Ifail,t0)
       IF (Ifail.LT.0) THEN
          Ilflag = 0
          lsfail = 1
          RETURN
       ENDIF
       rr = dsqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
       tet(1) = DACOS(x2(3)/rr)
       phi(1) = DATAN2(x2(2),x2(1))
!C
!C et on tourne -> on se decale sur la surface en phi et on cherche teta
!C pour avoir leI0 et B0 constants
!C
!C and one turns - > one shifts on surface in phi and one seeks teta to have
!C leI0 and B0 constant
!C
!      print *,'*** 8: find thetas around N. cap with I and Bm constan'
       dsreb = -dsreb

!       print *,'1,phi(1) =',1,phi(1)
!       print *,'1,tet(1)=',1,tet(1)/rad

       DO I = 2,Nder
        phi(I) = phi(I-1)+2.0*pi/Nder

!        print *,'i,phi(i) =',i,phi(i)

        Iflag_I = 0

!c       bkress -- I do not understand original code below (commented out)
!c       assumes Ilflag somewhere initialized to zero??
!c       initialize/set tetl
!*       IF (Ilflag.EQ.0) THEN
         tetl = tet(I-1)
         IF (I.GT.2) tetl = 2.0*tet(I-1)-tet(I-2)
         tet1 = tetl
!*       ELSE
!*        tetl = tet(I)
!*        tet1 = tetl
!*       ENDIF

!c       write(6,*)tetl
!c       read(5,*)
        leI1 = baddata
!C
107     CONTINUE
!c       initial point on earth's surface
        x1(1) = DSIN(tetl)*DCOS(phi(I))
        x1(2) = DSIN(tetl)*DSIN(phi(I))
        x1(3) = DCOS(tetl)
        Iflag = 0
        leI = baddata
!C
!c ********* bkress ??
            CALL CHAMP(x1,B,Bprev,Ifail,t0)
            IF (Ifail.LT.0) THEN
               Ilflag = 0
               lsfail = 1
               RETURN
            ENDIF
!        print *
!        print *,'i,tetl,phi(I),Bprev=',i,tetl/rad,phi(i)/rad,Bprev
!**********************
!C
!c       trace field line & get I
        DO J = 1,Nrebmax
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
            lsfail = 1
            RETURN
         ENDIF
         rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
!c        when we pass the northern mirror point,
!c        start calculating/integrating I
         IF (Bl.LT.B0) THEN
!c         initialize I (first time only, Iflag initially 0, set to 1 below)
          IF (Iflag .EQ. 0) THEN
            CALL CHAMP(x1,B,B1,Ifail,t0)
            IF (Ifail.LT.0) THEN
               Ilflag = 0
               lsfail = 1
               RETURN
            ENDIF
            leI = 0.5*DSQRT(1.0-Bl/B0)*(1.0+(Bl-B0)/(Bl-B1))
            Iflag = 1
          ELSE
            leI = leI+DSQRT(1.0-Bl/B0)
          ENDIF
         ENDIF
!c        break after we pass southern mirror point
         IF (Bl.GT.B0 .AND. Iflag.EQ.1) GOTO 103
         IF (rr2.LT.1.0) GOTO 103
         x1(1) = x2(1)
         x1(2) = x2(2)
         x1(3) = x2(3)
!c        write(6,*)J,Bl,B0,leI*ABS(dsreb)
        ENDDO
103     CONTINUE

!c Pourquoi?
!c       bkress -- commented out the following
!c       stepped inside of r=1 to find B > B0 (to get to mirror point)
        IF (rr2.LT.1.0) THEN
         leI = baddata
        ENDIF

!c       bkress -- modified the following
        IF (J.LT.Nrebmax .AND. rr2.GE.1.0) THEN
!*        IF (J.LT.Nrebmax) THEN
            CALL CHAMP(x1,B,B1,Ifail,t0)
            IF (Ifail.LT.0) THEN
               Ilflag = 0
               lsfail = 1
               RETURN
            ENDIF
            leI = leI+0.5*DSQRT(1.0-B1/B0)*(B0-Bl)/(Bl-B1)
            leI = leI*DABS(dsreb)
        ENDIF

!c       for first tetl only
!c       (Iflag_I set to 1 below)
        IF (Iflag_I .EQ.0) THEN
!c        ?? if we did not find southern mirror pt. (B0) ??
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

!c       if leI0 is between leI and leI1,
!c       we are done looking for correct theta
!c       (note, leI1 set = leI above,
!c       we do not goto 108 on the first time through)
        IF ((leI-leI0)*(leI1-leI0) .LT. 0.0) GOTO 108

        leI1 = leI
        tet1 = tetl
!c       *** set new tetl here ***
        IF (leI.LT.leI0) THEN
!c        if leI too small, move poleward
         tetl = tetl-dtet
        ElSE
!c        if leI too big, move equatorward
         tetl = tetl+dtet
        ENDIF

        IF (tetl.GT.pi .OR. tetl.LT.0.0) GOTO 108
        GOTO 107
108     CONTINUE

        tet(I) = 0.5*(tetl+tet1)

!c       done getting correct theta **********************
!*        print *,'i,tet(I)=',i,tet(I)/rad

!c	read(5,*)
        IF (J.GE.Nrebmax .AND. leI.GT.0.0) THEN
         Ilflag = 0
         lsfail = 1
         print *,'J.GE.Nrebmax .AND. leI.GT.0.0'
         RETURN
        ENDIF
!C
        x1(1) = DSIN(tet(I))*DCOS(phi(I))
        x1(2) = DSIN(tet(I))*DSIN(phi(I))
        x1(3) = DCOS(tet(I))

	CALL CHAMP(x1,B,Bl,Ifail,t0)
	IF (Ifail.LT.0) THEN
	   Ilflag = 0
           lsfail = 1
	   print *,'Ifail.LT.0'
	   RETURN
	ENDIF
	
*	print *,'Bl=',Bl
	
c       because of the method choosen (i.e. look for correct theta as we 
c       go around the polar cap at 1 RE) we have to restrict to L-shells 
c       that stay above 1RE. Here we check to see if the mirror point goes \
c       inside of 1RE.
	IF (Bl.LT.B0) THEN
	 Ilflag = 0
         lsfail = 1
	 print *,'Bl.LT.B0',Bl,B0
	 RETURN
	ENDIF
       ENDDO
C
C calcul de somme de BdS sur la calotte nord
C calculation of nap of BdS on the northern cap
C
***      print *,'*** 9: calc. sum of BdS on the northern cap'

       x1(1) = 0.0
       x1(2) = 0.0
       x1(3) = 1.0
       CALL CHAMP(x1,B,Bl,Ifail,t0)
       IF (Ifail.LT.0)THEN
	   Ilflag = 0
           lsfail = 1
	   RETURN
       ENDIF
       somme = Bl*pi*dtet*dtet/4.0
*       somme = pi*dtet*dtet/4.0
       DO I = 1,Nder
         tetl = 0.0
         DO J = 1,Ntet
	  tetl = tetl+dtet
	  IF (tetl .GT. tet(I)) GOTO 111
          x1(1) = DSIN(tetl)*DCOS(phi(I))
          x1(2) = DSIN(tetl)*DSIN(phi(I))
          x1(3) = DCOS(tetl)
          CALL CHAMP(x1,B,Bl,Ifail,t0)
          IF (Ifail.LT.0)THEN
	      Ilflag = 0
              lsfail = 1
	      RETURN
          ENDIF
c         here is error -- bkress
*          somme = somme+Bl*SIN(tetl)*dtet*2.0*pi/Nder ! Phi and not Lstar
          stetl=dsin(tetl)
          stetl2=stetl*stetl
          somme=somme-(b(1)*stetl2*dcos(phi(i))+b(2)*
     &       stetl2*dsin(phi(i))+b(3)*stetl*dcos(tetl))*dtet*2.0*pi/Nder
	 ENDDO
111      CONTINUE
       ENDDO

       Lstar = 2.0*pi*Bo/somme
       IF (Lm.LT.0.0) Lstar = -Lstar
       Ilflag = 1
C
*       print *,'somme,Lstar=',somme,Lstar

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
      REAL*8      y(3),Bxsm(3),Bl,t

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
         print *,'hit inner or outer boundary in L-shell calc.',status
         print *,'x,y,z,r=',
     &   y(1),y(2),y(3),dsqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
******         stop
         Ifail=-1
         return
      endif

c     b in sm coordinate system
      Bxsm(1)=bx/ffactor/ntg
      Bxsm(2)=by/ffactor/ntg
      Bxsm(3)=bz/ffactor/ntg
      Bl=dsqrt(bxsm(1)*bXsm(1)+bxsm(2)*bxsm(2)+bxsm(3)*bxsm(3))

c      if (sys.eq.2) print *,'need to rotate B to geo coordinates'

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
        REAL*8 xx(3),x2(3)
        REAL*8 B(3),Bl
        REAL*8 h
        REAL*8 xwrk(4,3)
        REAL*8 t0
C
C-----------------------------------------------------------------------
C
c        write(6,*)'sksyst'
c        write(6,*)xx(1),xx(2),xx(3),h
        CALL CHAMP(xx,B,Bl,Ifail,t0)
c	write(6,*)xx,B,Bl,Ifail
	IF (Ifail.LT.0) RETURN
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
	IF (Ifail.LT.0) RETURN
        xwrk(2,1) = h*B(1)/Bl
        xwrk(2,2) = h*B(2)/Bl
        xwrk(2,3) = h*B(3)/Bl
        x2(1) = xx(1)+xwrk(2,1)/2.0
        x2(2) = xx(2)+xwrk(2,2)/2.0
        x2(3) = xx(3)+xwrk(2,3)/2.0
c        write(6,*)x2(1),x2(2),x2(3),Bl
C
        CALL CHAMP(x2,B,Bl,Ifail,t0)
	IF (Ifail.LT.0) RETURN
        xwrk(3,1) = h*B(1)/Bl
        xwrk(3,2) = h*B(2)/Bl
        xwrk(3,3) = h*B(3)/Bl
        x2(1) = xx(1)+xwrk(3,1)
        x2(2) = xx(2)+xwrk(3,2)
        x2(3) = xx(3)+xwrk(3,3)
c        write(6,*)x2(1),x2(2),x2(3),Bl
C
        CALL CHAMP(x2,B,Bl,Ifail,t0)
	IF (Ifail.LT.0) RETURN
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
	IF (Ifail.LT.0) RETURN
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
c Modified S./ Bourdarie (July 2004)
c
C Routine to find mirror point of a trapped particle
C
       SUBROUTINE find_bmir(xsm,t0,alpha,Bposit,Bmir,xmin,ifail)
C
       IMPLICIT NONE
C
       integer  Nreb
       PARAMETER (Nreb = 50)
C
       integer  Ifail
       integer  Nrebmax
       real*8     rr
       real*8     xsm(3),t0,x1(3),x2(3)
       real*8     xmin(3)
       real*8     lati,longi,alti
       real*8     B(3),Bl,B0,B1,B3
       real*8     dsreb

       integer  I,J,ind
       real*8     Lb
       real*8     leI0
C
       real*8     pi,rad
       real*8     tt
       real*8     cste
C
       real*8     Bposit,Bmir
       real*8     sn2,sn2max,alpha
       
       REAL*8     baddata
       PARAMETER  (baddata=-1.d31)
C
*       COMMON /magmod/k_ext,k_l,kint
C
C
       pi = 4.0*ATAN(1.0)
       rad = pi/180.0
C
       Nrebmax = 20*Nreb
C
       leI0 = 0.0
C
       rr = DSQRT(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))
       tt = DACOS(xsm(3)/rr)
       Lb  = rr/DSIN(tt)/DSIN(tt)
c       write(6,*)'L bete ',Lb
C
       CALL CHAMP(xsm,B,B0,Ifail,t0)
       IF (Ifail.LT.0) THEN
          Bposit=baddata
	  Bmir=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
	  RETURN
       ENDIF
       Bmir = B0
       Bposit=B0
       sn2=DSIN(alpha)*DSIN(alpha)
       Sn2max=sn2
       cste=B0/sn2

c      set field line tracing step size
       dsreb = Lb/(Nreb*1.0)

C
C calcul du sens du depart 
C (calculation of the direction of the departure)

       CALL sksyst(-dsreb,xsm,x1,Bl,Ifail,t0)
       IF (Ifail.LT.0) THEN
	  Bmir=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
	  RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xsm,x2,Bl,Ifail,t0)
       IF (Ifail.LT.0) THEN
	  Bmir=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
	  RETURN
       ENDIF
       B3 = Bl

C
c      bkress -- modified to trace in direction specified by PA
c      rather than increasing field strength so that the mirror point location
c      returned is in the direction of the initial mirror point 
c      (bmir is the same either way).
*       IF (B1.GT.B3) THEN
*         dsreb = -dsreb
*       ENDIF
       IF (alpha.gt.(90.*rad)) THEN
         dsreb = -dsreb
       ENDIF

C
C calcul de la ligne de champ
C (calculation of the line of field)

       ind=0
       DO I = 1,3
         x1(I)  = xsm(I)
       ENDDO
C
c       write(6,*)dsreb
       xmin(1)=x1(1)
       xmin(2)=x1(2)
       xmin(3)=x1(3)
c      trace field line in direction of increasing magnetic field 
c      and record B at maximum B/B_mirror
       DO J = 1,Nrebmax
c        xx unchanged by sksyst, x2 gets new/advanced position
         CALL sksyst(dsreb,x1,x2,Bl,Ifail,t0)
         IF (Ifail.LT.0) THEN
	    Bmir=baddata
            xmin(1) = baddata
            xmin(2) = baddata
            xmin(3) = baddata
	    RETURN
         ENDIF
c        this is current Bl/B_mirror
	 sn2=Bl/cste
c        exit loop when B/B_mirror > 1
c        x1 is position immediately preceding
	 if (sn2 .GT.1.0) GOTO 20
         IF (sn2.GT.Sn2max) THEN
           xmin(1) = x2(1)
           xmin(2) = x2(2)
           xmin(3) = x2(3)
           Sn2max = sn2
	   Bmir=Bl
         ENDIF
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
20     CONTINUE
C
*       print *
*       print *,'position preceding mirror point (using dsreb steps)'
*       print *,'xmin=',xmin
*       print *,'|min|=',
*     & sqrt(xmin(1)*xmin(1)+xmin(2)*xmin(2)+xmin(3)*xmin(3))

       IF (J.GE.Nrebmax) THEN !open field line
	  Bmir=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
	  RETURN
       ENDIF
C           
C calcul de Bmirror point
C
c      now approach mirror point with smaller steps
       DO i=1,100
c        double to single
*         CALL sksyst(dsreb/100.0,xmin,x1,Bl,Ifail,t0)
         CALL sksyst(dsreb/100.0,xmin,x1,Bl,Ifail,t0)
         IF (Ifail.LT.0) THEN
	    Bmir=baddata
            xmin(1) = baddata
            xmin(2) = baddata
            xmin(3) = baddata
	    RETURN
         ENDIF
	 sn2=Bl/cste
c        again stop with x1 immediately preceding B/B_mirror > 1 location
	 if (sn2 .GT.1.0) GOTO 30
         xmin(1) = x1(1)
         xmin(2) = x1(2)
         xmin(3) = x1(3)
         Sn2max = sn2
	 Bmir=Bl
       ENDDO
30     CONTINUE

*       print *
*       print *,'final position immediately preceding mirror point 
*     & (using dsreb/100.0 steps)'
*       print *,'xmin=',xmin
*       print *,'|min|=',
*     & sqrt(xmin(1)*xmin(1)+xmin(2)*xmin(2)+xmin(3)*xmin(3))

       IF (xmin(1)*xmin(1)+xmin(2)*xmin(2)+
     &      xmin(3)*xmin(3).LT.1.0) THEN
	  Bmir=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
       ENDIF
C
100    CONTINUE
C
       END

C
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
       SUBROUTINE find_bmin(xsm,t0,Bloc,Bmin,xmin,ifail)
C
       IMPLICIT NONE
C
       integer  Nreb
       PARAMETER (Nreb = 50)
C
       integer  Ifail
       integer  Nrebmax
       real*8     rr
       real*8     xsm(3),t0,x0(3),x1(3),x2(3)
       real*8     xmin(3)
       real*8     B(3),Bl,B0,B1,B3
       real*8     dsreb
       integer  I,J,ind
       real*8     Lb
       real*8     pi,rad
       real*8     tt
       real*8     Bloc,Bmin
       REAL*8     aa,bb,smin
       REAL*8     baddata
       PARAMETER  (baddata=-1.d31)
C
       pi = 4.0*ATAN(1.0)
       rad = pi/180.0
C
       Nrebmax = 20*Nreb
C
       rr = DSQRT(xsm(1)*xsm(1)+xsm(2)*xsm(2)+xsm(3)*xsm(3))
       tt = DACOS(xsm(3)/rr)
       Lb  = rr/DSIN(tt)/DSIN(tt)
C
       CALL CHAMP(xsm,B,B0,Ifail,t0)
       IF (Ifail.LT.0) THEN
          Bloc=baddata
	  Bmin=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
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
       IF (Ifail.LT.0) THEN
	  Bmin=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
	  RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xsm,x2,Bl,Ifail,t0)
       IF (Ifail.LT.0) THEN
	  Bmin=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
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
         if (dabs(-0.5*bb).gt.dabs(aa)) then
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
         IF (Ifail.LT.0) THEN
            Bmin = baddata
            xmin(1) = baddata
            xmin(2) = baddata
            xmin(3) = baddata
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
         IF (Bl.GT.B0) GOTO 20
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
20     CONTINUE
       IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.0) THEN
        print *,'conjugate point inside 1RE surface'
       ENDIF
       IF (J.GE.Nrebmax) THEN !open field line
	  Bmin=baddata
          xmin(1) = baddata
          xmin(2) = baddata
          xmin(3) = baddata
	  RETURN
       ENDIF

c      now approach & hopfully pass bmin with smaller steps starting at x0
       DO i=1,200
         CALL sksyst(dsreb/100.0,x0,x1,Bl,Ifail,t0)
         IF (Ifail.LT.0) THEN
	    Bmin=baddata
            xmin(1) = baddata
            xmin(2) = baddata
            xmin(3) = baddata
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


