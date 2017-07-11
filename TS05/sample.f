
      program sample
	
c ************************************************************************* 
c
c Puts E&B fields on a cartisian grid. The output is an HDF file 
c for use with rbelt code. Note that the fields are scaled for use 
c with the rbelt code in here. The HDF output files are set up as follows:
c
c sds sets
c 1:(nx)x
c 2:(ny)y
c 3:(nz)z
c 4:(nt)t
c 5:(nt x nx x ny x nz)bx 
c 6:(nt x nx x ny x nz)by 
c 7:(nt x nx x ny x nz)bz 
c 8:(nt x nx x ny x nz)ex 
c 9:(nt x nx x ny x nz)ey 
c 10:(nt x nx x ny x nz)ez 
c	
c B.Kress
c
c *************************************************************************

      include 'hdf.inc'

      implicit none
c     declare parameter types
      integer nx,ny,nz,nt
      real xmin,xmax
      real ymin,ymax
      real zmin,zmax
      character* (*) basename
      real charge_sign
      real mass,dipole

c *************************************************************************
c BELOW ARE THE PARAMETERS TO SET:
c
c DEFINE RBELT GRID PARAMETERS:
c Use parameters below to define cartisian grid for use with rbelt code.

      parameter (nx=141, ny=141, nz=141, nt=1) !nt is #of time steps in output file
      parameter (xmin=-14., xmax=14.)
      parameter (ymin=-14., ymax=14.)
      parameter (zmin=-14., zmax=14.)

      parameter (basename='Tsy_Nov01_1-1')  !output file name

      parameter (charge_sign=-1)  !reverse charge_sign to run backwards in time
      parameter (mass=1.6726485e-24)  !proton mass
*      parameter (mass=9.109534e-28)  !electron mass
      parameter (dipole=30000.0)  !bzero dipole field (nT)
*      parameter (dipole=0.0)

c *************************************************************************

c		define output variables:
      real x(nx),y(ny),z(nz),t(nt)
      real bx(nx,ny,nz,nt),by(nx,ny,nz,nt),bz(nx,ny,nz,nt)
      real ex(nx,ny,nz,nt),ey(nx,ny,nz,nt),ez(nx,ny,nz,nt)

      COMMON /fields/bx,by,bz,ex,ey,ez
     
c		HDF stuff:
      integer sfstart, sfcreate
      integer sfwdata, sfendacc, sfend
      integer*4 file_id, sds_id, rank, status
      integer dim_sizes4d(4), start4d(4), edge4d(4),stride4d(4)
      integer dim_sizes1d(1), start1d(1), edge1d(1),stride1d(1)

      integer sffinfo
      integer sfselect
      integer n_datasets,n_file_attrs

c	...and a few more:
      integer i,j,k
      real nx_tmp,ny_tmp,nz_tmp !initial time
      real bx_tmp,by_tmp,bz_tmp
      real r2,r,r5
      real Re, q, c, m0, b0
      real tscal,bscal,escal,ffactor
      character*80 ofile
      character*4 string_out
      integer string_begin,string_end
      logical first

c		Define constants
c		to make sample more general we should scale the fields in rbelt code.
      Re = 6378.137e5
      q = 4.8032424e-10
      c = 2.99792458e10
      m0 = mass
      ffactor=charge_sign*q*Re/(c*c*m0)
      bscal=ffactor*1.e-5
      escal=ffactor*1.e6/c
      tscal = c/Re
      b0 = dipole

c		set up points
      nx_tmp = nx
      ny_tmp = ny
      nz_tmp = nz
      if (nx.eq.1) nx_tmp=2
      if (ny.eq.1) ny_tmp=2
      if (nz.eq.1) nz_tmp=2
      do i=1,nx
         x(i) = xmin + ((xmax-xmin)/(nx_tmp-1))*(i-1)
*         print *,'x(',i,')=',x(i)
      enddo
      do i=1,ny
         y(i) = ymin + ((ymax-ymin)/(ny_tmp-1))*(i-1)
*         print *,'y(',i,')=',y(i)
      enddo
      do i=1,nz
         z(i) = zmin + ((zmax-zmin)/(nz_tmp-1))*(i-1)
*         print *,'z(',i,')=',z(i)
      enddo

      first = .true.
      t(1)=0.0*tscal
      do k=1,nz
         do j=1,ny
            do i=1,nx
             r2=x(i)*x(i)+y(j)*y(j)+z(k)*z(k)
             r=sqrt(r2)
             r5=r2*r2*r
             call tsy_igrf(x(i),y(j),z(k),bx_tmp,by_tmp,bz_tmp,first)
             bx(i,j,k,1)=(bx_tmp+3.*b0*x(i)*z(k)/r5)*bscal
             by(i,j,k,1)=(by_tmp+3.*b0*y(j)*z(k)/r5)*bscal
             bz(i,j,k,1)=(bz_tmp+3.*b0*z(k)*z(k)/r5 - b0/r2/r)*bscal
             ex(i,j,k,1) = 0.
             ey(i,j,k,1) = 0.
             ez(i,j,k,1) = 0.
            enddo
         enddo
      enddo

c     then write out to HDF 
c     send 	x,y,z,t,bx,by,bz,ex,ey,ez out to HDF file ************

      rank=1
      dim_sizes1d(1) = nx
      start1d(1)=0
      stride1d(1)=1
      edge1d(1)=dim_sizes1d(1)

c     name and open the output file

      ofile=basename//'.hdf'
      file_id = sfstart(ofile,4)

      print *
      print *,'writing out to HDF file ',ofile

      sds_id=sfcreate(file_id,'x',DFNT_FLOAT32,rank,dim_sizes1d)
      status = sfwdata(sds_id,start1d,stride1d,edge1d,x)
      status = sfendacc(sds_id)
            
      dim_sizes1d(1) = ny
      edge1d(1)=dim_sizes1d(1)
            
      sds_id=sfcreate(file_id,'y',DFNT_FLOAT32,rank,dim_sizes1d)
      status = sfwdata(sds_id,start1d,stride1d,edge1d,y)
      status = sfendacc(sds_id)
            
      dim_sizes1d(1) = nz
      edge1d(1)=dim_sizes1d(1)
            
      sds_id=sfcreate(file_id,'z',DFNT_FLOAT32,rank,dim_sizes1d)
      status = sfwdata(sds_id,start1d,stride1d,edge1d,z)
      status = sfendacc(sds_id)
            
      dim_sizes1d(1) = 1
      edge1d(1)=dim_sizes1d(1)
            
      sds_id=sfcreate(file_id,'t',DFNT_FLOAT32,rank,dim_sizes1d)
      status = sfwdata(sds_id,start1d,stride1d,edge1d,t)
      status = sfendacc(sds_id)
            
c     output fields

      rank=4
      dim_sizes4d(1)=nx
      dim_sizes4d(2)=ny
      dim_sizes4d(3)=nz
      dim_sizes4d(4)=1
      do i=1,rank 
         stride4d(i)=1
         start4d(i)=0
         edge4d(i)=dim_sizes4d(i)
      enddo

      sds_id=sfcreate(file_id,'bx',DFNT_FLOAT32,rank,dim_sizes4d)
      status = sfwdata(sds_id,start4d,stride4d,edge4d,bx)
      status = sfendacc(sds_id)
            
      sds_id=sfcreate(file_id,'by',DFNT_FLOAT32,rank,dim_sizes4d)
      status = sfwdata(sds_id,start4d,stride4d,edge4d,by)
      status = sfendacc(sds_id)
            
      sds_id=sfcreate(file_id,'bz',DFNT_FLOAT32,rank,dim_sizes4d)
      status = sfwdata(sds_id,start4d,stride4d,edge4d,bz)
      status = sfendacc(sds_id)
            
      sds_id=sfcreate(file_id,'ex',DFNT_FLOAT32,rank,dim_sizes4d)
      status = sfwdata(sds_id,start4d,stride4d,edge4d,ex)
      status = sfendacc(sds_id)
            
      sds_id=sfcreate(file_id,'ey',DFNT_FLOAT32,rank,dim_sizes4d)
      status = sfwdata(sds_id,start4d,stride4d,edge4d,ey)
      status = sfendacc(sds_id)
            
      sds_id=sfcreate(file_id,'ez',DFNT_FLOAT32,rank,dim_sizes4d)
      status = sfwdata(sds_id,start4d,stride4d,edge4d,ez)
      status = sfendacc(sds_id)
      status = sfend(file_id)
            
      print *
      do i=1,nx
          print *,'x,bz=',x(i),bz(i,(ny+1)/2,(nz+1)/2,1)/bscal
      enddo

      end
	


	
