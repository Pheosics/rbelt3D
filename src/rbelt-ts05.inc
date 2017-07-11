c rbelt TS05 include file

      integer iopt
      real parmod(11),tilt0
      common /ts05params/parmod,tilt0,iopt

      character*(*) tsy_input
      logical tsy
      integer max_table
      parameter (tsy=.true.) !add Tsyganenko to IGRF
      parameter (tsy_input='TS05_input.dat') ! Tsyganenko/geopack input file name
      parameter (max_table=15000) ! must be > number of lines in T05_input file.
      
      integer num_table,year_tbl(max_table),doy_tbl(max_table),
     &hour_tbl(max_table),min_tbl(max_table),stsec_tbl,ndsec_tbl
      real pdyn(max_table),dst(max_table),bximf(max_table),
     &byimf(max_table),bzimf(max_table),tilt(max_table),
     &w1(max_table),w2(max_table),w3(max_table),
     &w4(max_table),w5(max_table),w6(max_table)
      common/ts05table/num_table,year_tbl,doy_tbl,hour_tbl,min_tbl,pdyn,
     &dst,bximf,byimf,bzimf,tilt,w1,w2,w3,w4,w5,w6,stsec_tbl,ndsec_tbl

