
c     rbelt-dst.inc include file

c     maximum number of lines (days) in dst data file
      integer maxdays,numdays
      parameter(maxdays=4)
      real*8 dummy(24),sfr(24*maxdays),dst(24*maxdays), tsave
      real*8 delt
      parameter(delt=3600*tfactor)
      common /dst/numdays,sfr,dummy,dst,tsave
