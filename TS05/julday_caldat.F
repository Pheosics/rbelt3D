
      FUNCTION julday(mm,id,iyyy)
      INTEGER julday,id,iyyy,mm,IGREG
      PARAMETER (IGREG=15+31*(10+12*1582)) 
c Gregorian Calendar adopted Oct. 15, 1582.
c In this routine julday returns the Julian Day Number that begins at noon ofthe calendar
c date specifed by month mm, day id, and year iyyy, all integer variables. Positive year
c signifes A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D.
      INTEGER ja,jm,jy
      jy=iyyy
      if (jy.eq.0) pause 'julday: there is no year zero'
      if (jy.lt.0) jy=jy+1
      if (mm.gt.2) then
         jm=mm+1
      else
         jy=jy-1
         jm=mm+13
      endif
      julday=365*jy+int(0.25d0*jy+2000.d0)+int(30.6001d0*jm)+id+1718995
      if (id+31*(mm+12*iyyy).ge.IGREG) then !Test whether to change to Gregorian Calendar.
         ja=int(0.01d0*jy)
         julday=julday+2-ja+int(0.25d0*ja)
      endif
      return
      END

      SUBROUTINE caldat(julian,mm,id,iyyy)
      INTEGER id,iyyy,julian,mm,IGREG
      PARAMETER (IGREG=2299161)
c Inverse of the function julday given above. Here julian is input as a Julian Day Number,
c and the routine outputs mm,id, and iyyy as the month, day, and year on which the specied
c Julian Day started at noon.
      INTEGER ja,jalpha,jb,jc,jd,je
      if(julian.ge.IGREG)then !Cross-over to Gregorian Calendar produces this correction.
         jalpha=int(((julian-1867216)-0.25d0)/36524.25d0)
         ja=julian+1+jalpha-int(0.25d0*jalpha)
      else if(julian.lt.0)then 
c        Make day number positive by adding integer number of Julian centuries,
c        then subtract them off at the end.
         ja=julian+36525*(1-julian/36525)
      else
         ja=julian
      endif
      jb=ja+1524
      jc=int(6680.0d0+((jb-2439870)-122.1d0)/365.25d0)
      jd=365*jc+int(0.25d0*jc)
      je=int((jb-jd)/30.6001d0)
      id=jb-jd-int(30.6001d0*je)
      mm=je-1
      if(mm.gt.12)mm=mm-12
      iyyy=jc-4715
      if(mm.gt.2)iyyy=iyyy-1
      if(iyyy.le.0)iyyy=iyyy-1
      if(julian.lt.0)iyyy=iyyy-100*(1-julian/36525)
      return
      END

