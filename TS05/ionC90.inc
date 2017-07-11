c ionC90.inc
	COMMON /ion/ psin(njp1,nk2p3),curn(njp1,nk2p3),
     $         pedern(njp1,nk2p3),
     $         halln(njp1,nk2p3),fluxn(njp1,nk2p3),e_avn(njp1,nk2p3),
     $         psis(njp1,nk2p3),curs(njp1,nk2p3),peders(njp1,nk2p3),
     $         halls(njp1,nk2p3),fluxs(njp1,nk2p3),e_avs(njp1,nk2p3)

	DIMENSION rionvar(nionvar)
	EQUIVALENCE (rionvar,psin)
