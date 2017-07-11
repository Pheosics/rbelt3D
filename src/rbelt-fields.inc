c    rbelt-fields.inc - particle E & B fields and spatial derivatives

      real b,bx,by,bz,ex,ey,ez,dbdx,dbdy,dbdz,dbxdx,dbxdy,dbxdz,
     &dbydx,dbydy,dbydz,dbzdx,dbzdy,dbzdz,dbdt,dbxdt,dbydt,dbzdt

      common /fields/b,bx,by,bz,ex,ey,ez,dbdx,dbdy,dbdz,
     &dbxdx,dbxdy,dbxdz,dbydx,dbydy,dbydz,dbzdx,dbzdy,dbzdz,dbdt,
     &dbxdt,dbydt,dbzdt

*      real vpotx,vpoty,vpotz
*      common /vectpot/vpotx,vpoty,vpotz

      real fnorm,b00
      common /fnormalize/fnorm,b00

      logical noe
      common /efield/noe

*      data bx /0./
*      data ex /0./
*      data ey /0./
*      data ez /0./
*      data bx /0./
*      data by /0./  
*      data bz=charge_sign*1.e-8
*      data b=1.e-8
*      data dbxdx /0./
*      data dbxdy /0./
*      data dbxdz /0./ 
*      data dbydx /0./
*      data dbydy /0./
*      data dbydz /0./
*      data dbzdx /0./ 
*      data dbzdy /0./
*      data dbzdz /0./
*      data dbdx /0./
*      data dbdy /0./
*      data dbdz /0./
*      data dbxdt /0./
*      data dbydt /0./
*      data dbzdt /0./
*      data dbdt /0./
