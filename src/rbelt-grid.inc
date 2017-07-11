c Set nx,ny,nz,nt, and nstep
c nx is x positions in Re in fields file
c ny is y positions in Re in fields file
c nz is z positions in Re in fields file
c nt is number of time steps in fields file
c nstep is number of time steps to read into rbelt x,y,z,t grid
c note, need > numx*numy*numz*nt*6*4*2 bytes of RAM to run this code.

      integer nx,ny,nz,nt,nstep

c     set grid parameters here
c     ****************************************************************

c     for a pure dipole or other analytically expressable fields
*      parameter (nx=1,ny=1,nz=1,nt=1)
*      parameter (nstep=1)

c     for time independent fields on a grid
*      parameter (nx=161,ny=161,nz=161,nt=1)
*      parameter (nstep=1)
*      parameter (nstep=3)

c     for time independent interpolated fields on a grid
c     (interpolates fields from time step grid positions 2 and 3
c     to 1, then only uses fields in 1)
*      parameter (nx=161,ny=161,nz=161,nt=1)
*      parameter (nstep=3)

c     for time dependent fields on a grid
      parameter (nx=241,ny=241,nz=121,nt=5)
      parameter (nstep=nt+2)

c     set grid boundaries here
      real*8 xgmin,ygmin,zgmin,xgmax,ygmax,zgmax
      parameter (xgmin=-6., xgmax=6.)
      parameter (ygmin=-6., ygmax=6.)
      parameter (zgmin=-3., zgmax=3.)

c     ****************************************************************

      real*8 xgr(nx),ygr(ny),zgr(nz),tgr(nstep)
      common /grid/xgr,ygr,zgr,tgr

      real*8 bxd(nx,ny,nz,nstep),byd(nx,ny,nz,nstep),
     1       bzd(nx,ny,nz,nstep),exd(nx,ny,nz,nstep),
     2       eyd(nx,ny,nz,nstep),ezd(nx,ny,nz,nstep)
*     3   lsg(nx,ny,nstep)

      integer nxy,nxyz,nxyzt
      parameter(nxy=nx*ny,nxyz=nx*ny*nz,nxyzt=nx*ny*nz*nstep)
      real*8 bxdv(nxyzt),bydv(nxyzt),bzdv(nxyzt),
     1       exdv(nxyzt),eydv(nxyzt),ezdv(nxyzt)

      equivalence (bxd(1,1,1,1),bxdv(1)),(byd(1,1,1,1),bydv(1)),
     1            (bzd(1,1,1,1),bzdv(1)),(exd(1,1,1,1),exdv(1)),
     2            (eyd(1,1,1,1),eydv(1)),(ezd(1,1,1,1),ezdv(1))

      real*8 dbxdxd(nx,ny,nz,nstep),dbxdyd(nx,ny,nz,nstep),
     &dbxdzd(nx,ny,nz,nstep),dbydxd(nx,ny,nz,nstep),
     &dbydyd(nx,ny,nz,nstep),dbydzd(nx,ny,nz,nstep),
     &dbzdxd(nx,ny,nz,nstep),dbzdyd(nx,ny,nz,nstep),
     &dbzdzd(nx,ny,nz,nstep),dbxdtd(nx,ny,nz,nstep),
     &dbydtd(nx,ny,nz,nstep),dbzdtd(nx,ny,nz,nstep)

      real*8 dbxdxv(nxyzt),dbxdyv(nxyzt),dbxdzv(nxyzt),
     &dbydxv(nxyzt),dbydyv(nxyzt),dbydzv(nxyzt),
     &dbzdxv(nxyzt),dbzdyv(nxyzt),dbzdzv(nxyzt),
     &dbxdtv(nxyzt),dbydtv(nxyzt),dbzdtv(nxyzt)

      equivalence (dbxdxd(1,1,1,1),dbxdxv(1)),
     1  (dbxdyd(1,1,1,1),dbxdyv(1)),
     4  (dbxdzd(1,1,1,1),dbxdzv(1)),(dbydxd(1,1,1,1),dbydxv(1)),
     5  (dbydyd(1,1,1,1),dbydyv(1)),(dbydzd(1,1,1,1),dbydzv(1)),
     6  (dbzdxd(1,1,1,1),dbzdxv(1)),(dbzdyd(1,1,1,1),dbzdyv(1)),
     7  (dbzdzd(1,1,1,1),dbzdzv(1)),(dbxdtd(1,1,1,1),dbxdtv(1)),
     8  (dbydtd(1,1,1,1),dbydtv(1)),(dbzdtd(1,1,1,1),dbzdtv(1))

      common /fgrid/bxd,byd,bzd,exd,eyd,ezd,
     &dbxdxd,dbxdyd,dbxdzd,dbydxd,dbydyd,dbydzd,dbzdxd,dbzdyd,
     &dbzdzd,dbxdtd,dbydtd,dbzdtd

*      common /fgrid/bxd,byd,bzd,exd,eyd,ezd

c     new arrays for derivatives ********************************
*****************************************************************
*      equivalence (dbxdxd(1,1,1,1),dbxdxv(1)),
*     &(dbxdyd(1,1,1,1),dbxdyv(1)),(dbxdzd(1,1,1,1),dbxdzv(1)),
*     &(dbydxd(1,1,1,1),dbydxv(1)),(dbydyd(1,1,1,1),dbydyv(1)),
*     &(dbydzd(1,1,1,1),dbydzv(1)),(dbzdxd(1,1,1,1),dbzdxv(1)),
*     &(dbzdyd(1,1,1,1),dbzdyv(1)),(dbzdzd(1,1,1,1),dbzdzv(1)),
*     &(dbxdtd(1,1,1,1),dbxdtv(1)),(dbydtd(1,1,1,1),dbydtv(1)),
*     &(dbzdtd(1,1,1,1),dbzdtv(1))

*****************************************************************

      real*8 dx,dy,dz
      common /metric/dx,dy,dz

      real*8 tzero1,tzero
      common /init_time/tzero1,tzero

c     sys = 1 for SM coordinates, 2 for GEO coordinates, ...
      integer sys
      common /coorsys/ sys
