c     rbelt-load_field.f - reads rbelt hdf files

************************************************************************

      subroutine load_field
     &(filename,filenum,firststep,laststep,gridstep)

c Reads E&B fields from an HDF file for use w/ rbelt code
c The HDF file SDS set format should be:
c
c 1:(nt)t
c 2:(nt x nx x ny x nz)bx 
c 3:(nt x nx x ny x nz)by 
c 4:(nt x nx x ny x nz)bz 
c 5:(nt x nx x ny x nz)ex 
c 6:(nt x nx x ny x nz)ey 
c 7:(nt x nx x ny x nz)ez 
c
c Opens filename & reads firststep up to & including laststep time steps
c (laststep-firststep+1 steps) of gridded fields with dims. nx x ny x nz.
c The fields are put into arrays bxd(nx,ny,nz,nstep), byd(nx,,ny,... etc.
c defined in grid.inc. The fields are read into bxd, byd, etc. starting 
c with location e.g. bxd(nx,ny,nz,gridstep) up to 
c bxd(nx,ny,nz,(gridstep+laststep-firststep)). filenum not used here.
c
c E.g., To read all fields from an HDF field file with nstep time steps
c put firststep=1, laststep=nstep, & gridstep=1.
c
c Note that: laststep must be > firststep & (gridstep+laststep-firststep)
c must be <= nstep. B.Kress -- 2/27/06

************************************************************************
 
      implicit none
      include 'hdf.inc' !needed for hdf commands
      include 'dffunc.inc' !defines hdf functions
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'

      integer i,j,k,l
      integer filenum,firststep,laststep,gridstep,index
      character*(*) filename
      integer lnblnk

      integer rank, secid, secindex
      integer fileid, status
      integer dim1d(1), edges1d(1), start1d(1), stride1d(1)
      integer dim4d(4), edges4d(4), start4d(4), stride4d(4)

      real tgrtemp(nstep),bxdvtemp(nxyzt),bydvtemp(nxyzt)
      real bzdvtemp(nxyzt),exdvtemp(nxyzt),eydvtemp(nxyzt)
      real ezdvtemp(nxyzt)

      real*8 bxsave(nx,ny,nz,2),bysave(nx,ny,nz,2),bzsave(nx,ny,nz,2)
      real*8 exsave(nx,ny,nz,2),eysave(nx,ny,nz,2),ezsave(nx,ny,nz,2)

!      print *
!      print *,'*** in subroutine hdffieldread ***'
!      print *,'nx,ny,nz,nstep=',nx,ny,nz,nstep
!      print *,'opening ',filename(1:lnblnk(filename))
      fileid = sfstart(filename(1:lnblnk(filename)),DFACC_READ)     ! sd_id 

*      dim1d(1) = nt
      start1d(1) = firststep-1 !first step in HDF file is step zero
      stride1d(1) = 1
      edges1d(1) = (laststep-firststep+1)
      secindex = sfn2index(fileid,'t')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start1d,stride1d,edges1d,tgrtemp(gridstep))
      if (status .eq. -1) then
         print*,'HDF read failed at tgr'
         stop
      endif
      status = sfendacc(secid)

!      print *,'read time steps:'
      i=gridstep
      if (i.ne.1) tgrtemp(i) = tgrtemp(i)! + 8640000
!      print *,'i,tgr(i)',i,tgrtemp(i)
      do i=gridstep+1,gridstep+(laststep-firststep+1)-1
*         print *,'i,tgr(i),delta',i,tgr(i),tgr(i)-tgr(i-1)
         tgrtemp(i) = tgrtemp(i)! + 8640000
!         print *,'i,tgr(i)',i,tgrtemp(i)
      enddo

c     Select and read in field data

      index=nxyz*(gridstep-1)+1

      rank=4             ! nx, ny, nz and nt are the 4
      dim4d(1)=nx        ! dim_sizes(1)
      dim4d(2)=ny        ! dim_sizes(2)
      dim4d(3)=nz        ! dim_sizes(3)

      do i=1,rank-1
         start4d(i)=0
         stride4d(i)=1
         edges4d(i)=dim4d(i)
      enddo
      start4d(4)=firststep-1 !first step in HDF file is step zero
      stride4d(4)=1
      edges4d(4)=(laststep-firststep+1)

      secindex = sfn2index(fileid,'bx')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start4d,stride4d,edges4d,bxdvtemp(index))
      if (status .eq. -1) then
         print*,'HDF read failed at bxd read'
         stop
      endif
      status = sfendacc(secid)

      secindex = sfn2index(fileid,'by')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start4d,stride4d,edges4d,bydvtemp(index))
      status = sfendacc(secid)

      secindex = sfn2index(fileid,'bz')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start4d,stride4d,edges4d,bzdvtemp(index))
      status = sfendacc(secid)

      secindex = sfn2index(fileid,'ex')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start4d,stride4d,edges4d,exdvtemp(index))
      status = sfendacc(secid) 

      secindex = sfn2index(fileid,'ey')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start4d,stride4d,edges4d,eydvtemp(index))
      status = sfendacc(secid)

      secindex = sfn2index(fileid,'ez')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start4d,stride4d,edges4d,ezdvtemp(index))
      status = sfendacc(secid)

c     close hdf file
      status = sfend(fileid)
      if (status .eq. -1) then
         print*,'HDF file close failed'
      endif

      if (filenum.ne.1) then
        tgrtemp(1) = tgr(1)
        tgrtemp(2) = tgr(2)
        do l=1,2
          do k=1,nz
            do j=1,ny
              do i=1,nx
                bxsave(i,j,k,l)=bxd(i,j,k,l)
                bysave(i,j,k,l)=byd(i,j,k,l)
                bzsave(i,j,k,l)=bzd(i,j,k,l)
                exsave(i,j,k,l)=exd(i,j,k,l)
                eysave(i,j,k,l)=eyd(i,j,k,l)
                ezsave(i,j,k,l)=ezd(i,j,k,l)  
              enddo
            enddo
          enddo
        enddo
      endif

      tgr = tgrtemp
      bxdv = bxdvtemp
      bydv = bydvtemp
      bzdv = bzdvtemp
      exdv = exdvtemp
      eydv = eydvtemp
      ezdv = ezdvtemp


      if (filenum.ne.1) then
        do l=1,2
          do k=1,nz
            do j=1,ny
              do i=1,nx
                bxd(i,j,k,l)=bxsave(i,j,k,l)
                byd(i,j,k,l)=bysave(i,j,k,l)
                bzd(i,j,k,l)=bzsave(i,j,k,l)
                exd(i,j,k,l)=exsave(i,j,k,l)
                eyd(i,j,k,l)=eysave(i,j,k,l)
                ezd(i,j,k,l)=ezsave(i,j,k,l)
              enddo
            enddo
          enddo
        enddo
      endif


*       print *
*       do j=1,nt
*          print *,filename,j
*          do i=1,nx
*            print *,xgr(i),1000.0*
*     &      eyd(i,(ny+1)/2,(nz+1)/2,j)/ffactor/vmsvcm
*          enddo
*       enddo

      return
      end

************************************************************************

      subroutine load_grid(filename)

      implicit none
      include 'hdf.inc'
      include 'dffunc.inc'
      include 'rbelt-grid.inc'

      integer i,j,k,l
      integer firststep,laststep,gridstep,index
      character*(*) filename
      character*80 ofile
      integer lnblnk

      integer rank, secid, secindex
      integer fileid, status
      integer dim1d(1), edges1d(1), start1d(1), stride1d(1)

      real xgrtemp(nx),ygrtemp(ny),zgrtemp(nz)

!      print *
!      print *,'*** in subroutine load_grid ***'
!      print *,'nx,ny,nz,nstep=',nx,ny,nz,nstep
!      print *,' opening ',filename(1:lnblnk(filename))
      fileid = sfstart(filename(1:lnblnk(filename)),DFACC_READ)

c     Select and read in grid data

*      rank = 1
      start1d(1) = 0
      stride1d(1) = 1

*      dim1d(1) = nx
      edges1d(1) = nx
      secindex = sfn2index(fileid,'x') ! sds_index
      secid = sfselect(fileid,secindex) ! sds_id
      status = sfrdata(secid,start1d,stride1d,edges1d,xgrtemp)
      if (status .eq. -1) then
         print*,'HDF read failed at xgr'
         stop
      endif
      status = sfendacc(secid)

*      dim1d(1) = ny
      edges1d(1) = ny
      secindex = sfn2index(fileid,'y')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start1d,stride1d,edges1d,ygrtemp)
      if (status .eq. -1) then
         print*,'HDF read failed at ygr'
         stop
      endif
      status = sfendacc(secid)

*      dim1d(1) = nz
      edges1d(1) = nz
      secindex = sfn2index(fileid,'z')
      secid = sfselect(fileid,secindex)
      status = sfrdata(secid,start1d,stride1d,edges1d,zgrtemp)
      if (status .eq. -1) then
         print*,'HDF read failed at zgr'
         stop
      endif
      status = sfendacc(secid)

c     close hdf file
      status = sfend(fileid)
      if (status .eq. -1) then
         print*,'HDF file close failed'
      endif
      
      xgr = xgrtemp
      ygr = ygrtemp
      zgr = zgrtemp

      return
      end

************************************************************************

      subroutine load_ut(basename)

c     read ut file
      implicit none
      include 'rbelt-ut.inc'
      character*80 filename,utfile,basename
      integer lnblnk
!      print *
!      print *,'*** in subroutine read_ut ***'
c     open ut ut file
      filename=utfile(basename)
!      print *,'opening ',filename(1:lnblnk(filename))
      open (12,file=filename(1:lnblnk(filename)))
      read(12,*)year0,doy0,hour0,min0,sec0
      close(12)
!      print *,'year0,doy0,hour0,min0,sec0=',year0,doy0,hour0,min0,sec0
      return
      end






