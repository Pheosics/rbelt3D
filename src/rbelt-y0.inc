c     rbelt-y0.inc - particle distribution state vector
c     set total number of particles here

      integer num_particles,num_pdata
      parameter(num_particles=1000000)
      parameter(num_pdata=10)
      real*8 y0(num_pdata,num_particles)
      integer int_y0(2,num_particles)
      common /dist/ y0,int_y0

