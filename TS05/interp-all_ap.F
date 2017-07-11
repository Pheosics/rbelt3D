	subroutine interp1d(datao,data1d,nl,nm,nn,nlnmnn,npts,nvar)
	include 'param.inc'
	parameter ( nbig = ni*njp2*nkp1)
c     in the following changed nlnmnn -> nlnmnn*(nvar+3) - bkress
*	real datao(nvar,npts),data1d(nlnmnn)
	real datao(nvar,npts),data1d(nlnmnn*(nvar+3))
	real data(ni,njp2,nkp2)
	real dataex(nbig)
	equivalence (data,dataex)
	include 'intco1d.inc'
c
	nlnm = nl*nm*(nvar+3)
	nv3 = nvar+3
	nv3nl = nv3*nl
*	write (6,*) (data1d(1+14*(i-1)),i=1,12)
c
        do 500 kvar=1,nvar
c
	do 40 k=1,nkp1
	do 40 j=1,njp2
	do 40 i=1,ni
	 data(i,j,k) = data1d(nlnm*(k-1)+nv3nl*(j-1)+(i-1)*nv3 +kvar)
*	data(i,j+1,k+1)  = float(i)*1.e-24
   40	continue

	do 42 j=1,njp2
	do 42 i = 1,ni
	  data(i,j,nkp2) = data(i,j,2)
   42	continue
c
	do 60 k=1,npts
c
	izero = icube(k)
	if (izero .gt. nlnmnn .or. izero .lt. 1) then
	   write(6,*) 'bad icube at ',i,j,k,izero
	endif
c
	datao(kvar,k) = 0.0
	datao(kvar,k) = dataex(izero)   * a(k,1)
	datao(kvar,k) = dataex(izero+ni) * a(k,2) 
     $			+ datao(kvar,k)
	datao(kvar,k) = dataex(izero+njp2*ni) * a(k,3) 
     $			+ datao(kvar,k)
	datao(kvar,k) = dataex(izero+njp3*ni) * a(k,4) 
     $			+ datao(kvar,k)
	datao(kvar,k) = dataex(izero+1) * a(k,5) 
     $			+ datao(kvar,k)
	datao(kvar,k) = dataex(izero+ni+1) * a(k,6) 
     $			+ datao(kvar,k)
	datao(kvar,k) = dataex(izero+njp2*ni+1) * a(k,7) 
     $			+ datao(kvar,k)
	datao(kvar,k) = dataex(izero+njp3*ni+1) * a(k,8) 
     $			+ datao(kvar,k)
c
c
   60	continue
 500    continue
cc
	return
	end










