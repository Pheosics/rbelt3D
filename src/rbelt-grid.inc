c Set nx,ny,nz,nt, and nt
c nx is number of x grid points
c ny is number of y grid points
c nz is number of z grid points
c nt is number of time grid points
c need > numx*numy*numz*nt*6*4*2 bytes of RAM to run this code.

      integer nx,ny,nz,nt

c     set grid parameters here
c     ****************************************************************

c     for a pure dipole or other analytically expressable fields
*      parameter (nx=1,ny=1,nz=1,nt=1)

c     for time independent fields on a grid
c     for single static field snapshot (one time grid position)
*      parameter (nx=241,ny=241,nz=241,nt=1)

c     for time independent fields on a grid
c     for static field snapshot in rbelt time grid position 3,
c     interpolated from fields in rbelt time grid positions 1 and 2.
*      parameter (nx=241,ny=241,nz=121,nt=7)

c     for time independent interpolated fields on a grid
c     (interpolates fields from time step grid positions 2 and 3
c     to 1, then only uses fields in 1)
*      parameter (nx=161,ny=161,nz=161,nt=1)

c     for time dependent fields on a grid
      parameter (nx=241,ny=241,nz=121,nt=7)

c     set grid boundaries here
      real xgmin,ygmin,zgmin,xgmax,ygmax,zgmax
      parameter (xgmin=-6., xgmax=6.)
      parameter (ygmin=-6., ygmax=6.)
      parameter (zgmin=-3., zgmax=3.)

c     ****************************************************************

      real xgr(nx),ygr(ny),zgr(nz),tgr(nt)
      common /grid/xgr,ygr,zgr,tgr

      real bxd(nx,ny,nz,nt),byd(nx,ny,nz,nt),
     1   bzd(nx,ny,nz,nt),exd(nx,ny,nz,nt),
     2   eyd(nx,ny,nz,nt),ezd(nx,ny,nz,nt)
*     3   lsg(nx,ny,nt)

      integer nxy,nxyz,nxyzt
      parameter(nxy=nx*ny,nxyz=nx*ny*nz,nxyzt=nx*ny*nz*nt)
      real bxdv(nxyzt),bydv(nxyzt),bzdv(nxyzt),
     1   exdv(nxyzt),eydv(nxyzt),ezdv(nxyzt)
      equivalence (bxd(1,1,1,1),bxdv(1)),(byd(1,1,1,1),bydv(1)),
     1   (bzd(1,1,1,1),bzdv(1)),(exd(1,1,1,1),exdv(1)),
     2   (eyd(1,1,1,1),eydv(1)),(ezd(1,1,1,1),ezdv(1))

*      common /fgrid/bxd,byd,bzd,exd,eyd,ezd

c     new arrays for derivatives ********************************
*****************************************************************
      real dbxdxd(nx,ny,nz,nt),dbxdyd(nx,ny,nz,nt),
     &dbxdzd(nx,ny,nz,nt),dbydxd(nx,ny,nz,nt),
     &dbydyd(nx,ny,nz,nt),dbydzd(nx,ny,nz,nt),
     &dbzdxd(nx,ny,nz,nt),dbzdyd(nx,ny,nz,nt),
     &dbzdzd(nx,ny,nz,nt),dbxdtd(nx,ny,nz,nt),
     &dbydtd(nx,ny,nz,nt),dbzdtd(nx,ny,nz,nt)

      real dbxdxv(nxyzt),dbxdyv(nxyzt),dbxdzv(nxyzt),
     &dbydxv(nxyzt),dbydyv(nxyzt),dbydzv(nxyzt),
     &dbzdxv(nxyzt),dbzdyv(nxyzt),dbzdzv(nxyzt),
     &dbxdtv(nxyzt),dbydtv(nxyzt),dbzdtv(nxyzt)

      equivalence (dbxdxd(1,1,1,1),dbxdxv(1)),
     &(dbxdyd(1,1,1,1),dbxdyv(1)),(dbxdzd(1,1,1,1),dbxdzv(1)),
     &(dbydxd(1,1,1,1),dbydxv(1)),(dbydyd(1,1,1,1),dbydyv(1)),
     &(dbydzd(1,1,1,1),dbydzv(1)),(dbzdxd(1,1,1,1),dbzdxv(1)),
     &(dbzdyd(1,1,1,1),dbzdyv(1)),(dbzdzd(1,1,1,1),dbzdzv(1)),
     &(dbxdtd(1,1,1,1),dbxdtv(1)),(dbydtd(1,1,1,1),dbydtv(1)),
     &(dbzdtd(1,1,1,1),dbzdtv(1))

      common /fgrid/bxd,byd,bzd,exd,eyd,ezd,dbxdxd,dbxdyd,dbxdzd,
     &dbydxd,dbydyd,dbydzd,dbzdxd,dbzdyd,dbzdzd,dbxdtd,dbydtd,dbzdtd

*****************************************************************

      real dx,dy,dz
      common /gdelta/dx,dy,dz

      real tzero1,tzero
      common /init_time/tzero1,tzero

c     sys = 1 for SM coordinates, 2 for GEO coordinates, ...
      integer sys
      common /coorsys/sys
