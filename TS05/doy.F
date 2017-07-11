************************************************************************

      subroutine doy_mmdd (mode, year, doy, mm, dd, status)

C Converts day-of year to months and days within month (mode = 1)
C Converts months and days to day-of-year (mode = -1)
C Valid for years 1901 - 2099 inclusive

C Exit status
C             0  normal exit
C             1  year is outside valid range
C             2  illegal value for mode
C             3  this part of code not checked


      integer*4       mode
      integer*4       year
      integer*4       doy
      integer*4       mm
      integer*4       dd
      integer*4       status
      integer*4       i

      integer*4       month(12) /31,28,31,30,31,30,31,31,30,31,30,31/

C Initializations

      status = 0

      if ((year .le. 1900) .or. (year .ge. 2100)) then
        status = 1
        go to 90
      end if

      month(2) = 28
      if (mod(year,4) .eq. 0) month(2) = 29


C Enter loops for calculation

      if (mode .eq. 1) then
        i = 1
        dd = doy
        do while (dd .gt. month(i))
          dd = dd - month(i)
          i = i + 1
        end do
        mm = i

      else if (mode. eq. -1) then

        doy = 0
        i = 1
        do while (i .lt. mm)
          doy = doy + month(i)
          i = i + 1
        end do
        doy = doy + dd

      else
        status = 2
        go to 90

      end if

C Exit

   90 continue

      return
      end

