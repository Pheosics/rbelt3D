************************************************************************

	SUBROUTINE INT2STR(CAOUT,IVALUE,IBEGIN,IEND)

	character*(*) caout
	integer iValue, iBegin, iEnd
	integer i,iStart,iLast
	character*20 caDum

	write(caDum,10) iValue
 10     format(i8)

	iStart = 1
	do while( caDum(iStart:iStart) .eq. ' ' )
	  iStart = iStart + 1
	enddo
	iLast = lnblnk(caDum)

	caOut(1:iLast-iStart+1) = caDum(iStart:iLast)
        iBegin = 1
	iEnd = iLast - iStart + 1

	return

	end

c this is an intrinsic
************************************************************************
*
*	integer function lnblnk(str)
*
*C
*C Description:
*C
*C This function returns the location of the last non-blank
*C character in a string. If all the characters are blank
*C the function returns 0.
*C
*C Important note: SUN FORTRAN provides this function. SUNS version is probably
*C more efficient than this one so use SUNs if you are using the SUN
*C compiler. This routine is provided only for those FORTRANs that 
*C do not have this function.
*C
*C User interface:
*C
*C Name         in/out type          Structure  Meaning
*C lnblnk       out    integer       scalar     The position of the last       
*C                                              non-blank character in the    
*C                                              input string. (Zero is returned
*C                                              if the entire string is blank)
*C str          in     character(*)  scalar     The string in which the last
*C                                              non-blank character is to be
*C                                              found.
*C
*C Errors:
*C
*C There should be no errors returned from this routine.
*C
*      implicit none
*      integer limit,i
*      logical done
*      character*(*) str
*      limit=1
*      i=len(str)
*      done=.false.
*      do while (.not. done)
*         if (str(i:i) .eq. ' ') then
*            i=i-1
*            if (i .lt. limit) then
*               done=.true.
*               lnblnk=0
*            end if
*         else
*            done=.true.
*            lnblnk=i
*         endif
*      enddo
*      return
*      end
*
************************************************************************

      character*80 function gridfile(basename,filenum)

      implicit none
      character*(*) basename
      integer filenum,lnblnk

      gridfile=basename(1:lnblnk(basename))//'grid.hdf'

      return
      end


************************************************************************

      character*80 function fieldfile(basename,filenum)

      implicit none
      character*(*) basename
      integer filenum,lnblnk,string_begin,string_end
      character*8 string_out

      call int2str(string_out,filenum,string_begin,string_end)
      fieldfile=basename(1:lnblnk(basename))//'field-'//
     &string_out(1:string_end)//'.hdf'

      return
      end

************************************************************************

      character*80 function distfile(basename,filenum)

      implicit none
      character*(*) basename
      integer filenum,lnblnk,string_begin,string_end
      character*8 string_out

      call int2str(string_out,filenum,string_begin,string_end)
      distfile=basename(1:lnblnk(basename))//'dist-'//
     &string_out(1:string_end)//'.dat'

      return
      end

************************************************************************

      character*80 function fluxfile(basename,filenum)

      implicit none
      character*(*) basename
      integer filenum,lnblnk,string_begin,string_end
      character*8 string_out

      call int2str(string_out,filenum,string_begin,string_end)
      fluxfile=basename(1:lnblnk(basename))//'flux-'//
     &string_out(1:string_end)//'.dat'

      return
      end

************************************************************************

      character*80 function prcpfile(basename,filenum)

      implicit none
      character*(*) basename
      integer filenum,lnblnk,string_begin,string_end
      character*8 string_out

      call int2str(string_out,filenum,string_begin,string_end)
      prcpfile=basename(1:lnblnk(basename))//'prcp-'//
     &string_out(1:string_end)//'.dat'

      return
      end
************************************************************************

      character*80 function conefile(basename,filenum)

      implicit none
      character*(*) basename
      integer filenum,lnblnk,string_begin,string_end
      character*8 string_out

      call int2str(string_out,filenum,string_begin,string_end)
      conefile=basename(1:lnblnk(basename))//'cone-'//
     &string_out(1:string_end)//'.dat'

      return
      end

************************************************************************

      character*80 function infofile(basename,filenum)

      implicit none
      character*(*) basename
      integer filenum,lnblnk,string_begin,string_end
      character*8 string_out

      call int2str(string_out,filenum,string_begin,string_end)
      infofile=basename(1:lnblnk(basename))//'info-'//
     &string_out(1:string_end)//'.txt'

      return
      end

************************************************************************

      character*80 function initfile(basename,filenum)

      implicit none
      character*(*) basename
      integer filenum,lnblnk,string_begin,string_end
      character*8 string_out

      print *,'*** in function initfile ***'

      call int2str(string_out,filenum,string_begin,string_end)
      initfile=basename(1:lnblnk(basename))//'init-'//
     &string_out(1:string_end)//'.dat'

      return
      end

************************************************************************

      character*80 function utfile(basename)

      implicit none
      character*(*) basename
      integer lnblnk
      utfile=basename(1:lnblnk(basename))//'ut.txt'
      return
      end

c*****************************************************************
c                    SUBROUTINE STR2INT                          *
c from:
c http://eosweb.larc.nasa.gov/PRODOCS/erbe/read_software/rds4gnshdf.f
c*****************************************************************
c Name - str2int                     Type - subroutine
c Version - 1.0   date - 07/28/93   Programmer- Nichele Brown(SAIC)
c Purpose - This subroutine converts a character string to an integer
c Input parameter - buffer : the character string to be converted
c Output parameter - value : the converted integer value
c                    iret : return flag
c                           0 - conversion was a success
c                           1 - no convertible string
c                           2 - string contained only spaces
c Key Local Parameters - subval : holds the value being converted
c                        mark : holds the position within the string
c                        nest : holds a single char. conversion
c                        loopcnt : loop counter
c                        spcnt : number of spaces in the string
c Subroutines called - none.
c*****************************************************************
      subroutine str2int(buffer, value, start, iret)
*      character buffer*80
      character buffer*(*)
      integer subval, value, sign
      integer mark, next, loopcnt, spcnt, iret, start
      spcnt = 0
      value = 0
      sign = 1
      iret = 0
c ****** validate that the string represents a numerical value ******
      do 100 loopcnt = start, len(buffer)
        if (buffer(loopcnt:loopcnt) .eq. ' ') spcnt = spcnt + 1
*        print *
*        print *,'buffer(loopcnt:loopcnt)=',buffer(loopcnt:loopcnt)
        if ((buffer(loopcnt:loopcnt) .lt. '0') .or.
     &  (buffer(loopcnt:loopcnt) .gt. '9')) then
          if (buffer(loopcnt:loopcnt) .ne. ' ') iret = 1
        endif
100   continue
      mark = start
      subval = 0
      if (spcnt .eq. len(buffer)) then
        iret = 2
      else
c ******  skip the leading spaces of the string ******
        do while (buffer(mark:mark) .eq. ' ')
          mark = mark + 1
        end do
c ****** convert the number from a string to a real ******
        if (buffer(mark:mark) .eq. '-') then
          sign = -1
          mark = mark + 1
        endif
        do while (buffer(mark:mark) .ne. ' ' .and.
     &            ichar(buffer(mark:mark)) .ne. 0)
          next = ichar(buffer(mark:mark)) - 48
          subval = subval * 10.0 + next
          mark = mark + 1
        end do
        do while (mark .le. len(buffer))
          if (buffer(mark:mark) .ne. ' ' .and. 
     &        ichar(buffer(mark:mark)) .ne. 0) iret = 1
          mark = mark + 1
        end do
      endif
      value = subval * sign
1100  return
      end

********************************************************************************
c although this is not related to string management, this is a convenient place for now 

c Numerical Rec. Portable rand# gen.

      real*8 FUNCTION ran0(idum) 
      INTEGER idum,IA,IM,IQ,IR,MASK 
*      REAL ran0,AM 
      REAL*8 AM 
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM, 
     &IQ=127773,IR=2836,MASK=123459876)  
      INTEGER k 

c Minimal  random number generator of Park and Miller. 
c Returns a uniform random deviate between 0.0 and 1.0. 
c Set or reset idum to any integer value (except the unlikely value MASK) 
c to initialize the sequence; idum must not be altered between calls for 
c successive deviates in a sequence. 

      idum=ieor(idum,MASK) 
      k=idum/IQ 
      idum=IA*(idum-k*IQ)-IR*k 
      if (idum.lt.0) idum=idum+IM 
      ran0=AM*idum 
      idum=ieor(idum,MASK)
      return 
      END





