      program main

      INTEGER hdftake, ierr
      INTEGER numx,numy,numz,npts
      LOGICAL interp
      PARAMETER (numx=101, numy=101, numz=101)
      PARAMETER (xmin=-10., xmax=10.)
      PARAMETER (ymin=-10., ymax=10.)
      PARAMETER (zmin=-10., zmax=10.)
      REAL xlim(6),allcart(14,numx,numy,numz)
      REAL points(3,numx,numy,numz) 
      REAL time
      CHARACTER(256) fname



      write(6,*) 'Starting HDF test function'
      fname='test.hdf'
      interp = .true.
      xlim(1)=xmin
      xlim(2)=xmax
      xlim(3)=ymin
      xlim(4)=ymax
      xlim(5)=zmin
      xlim(6)=zmax
      mstep=0
      npts=numx*numy*numz

      ierr = hdftake(fname,xlim,interp,allcart,points,npts,time)

      END
