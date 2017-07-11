	parameter (num1d=1030301)

	common/intrp/ a(num1d,8) , icube(num1d)
	common/carts/ 
     $		x2mid(ni,njp2),y2mid(ni,njp2),phimid(nkp2),
     $		x3mid(ni,njp2,nkp2),y3mid(ni,njp2,nkp2),
     $		z3mid(ni,njp2,nkp2),
     $	        psimid_in(njp2),psimid_out(njp2)
