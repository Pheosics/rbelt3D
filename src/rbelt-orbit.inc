c     rbelt-orbit.inc - rbelt-orbit.f include file

      integer max,num_tbl
      parameter(max=14100)

      integer year_tbl(max),doy_tbl(max),hour_tbl(max),min_tbl(max),
     &sec_tbl(max)
      real ftime,gcdist,lshell(max),mlat(max),mdl(max),alt(max),
     &lon(max),lat(max),gcr(max)
      common /table/num_tbl,doy_tbl,year_tbl,hour_tbl,min_tbl,sec_tbl,
     &ftime,gcdist,lshell,mlat,mdl,alt,lon,lat,gcr

      integer ent_index,ext_index,lwr_index,upr_index
      common /indcies/ent_index,ext_index,lwr_index,upr_index

c     for now we put flight start time here ...
      integer year0tbl,doy0tbl,hour0tbl,min0tbl,sec0tbl
      parameter(year0tbl=2003,doy0tbl=302)
      parameter(hour0tbl=21,min0tbl=0,sec0tbl=0)

      integer start_flag
      real pdiff
      common /stflag/pdiff,start_flag
      data start_flag/1/

      real max_mlat,min_mlat,energy,altitude
      integer strtyear,strtdoy,strthour,strtmin,strtsec
      common /inlist/max_mlat,min_mlat,energy,altitude,strtyear,strtdoy,
     &strthour,strtmin,strtsec
