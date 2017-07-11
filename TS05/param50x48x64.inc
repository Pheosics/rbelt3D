c param.inc
#include "lparam.inc"
#ifndef DOUBLEMESH
      PARAmeter (NI=50,NJ=48,NK=64)
#else
      PARAmeter (NI=100,NJ=48,NK=64)
#endif
      PARAmeter (
     Q    NIP1=NI+1, NIP2=NI+2, NIP3=NI+3, NIP4=NI+4,  
     Q    NIM1=NI-1, NIM2=NI-2, NIM3=NI-3, NIM4=NI-4)  
      PARAmeter (                        
     Q     NJP1=NJ+1, NJP2=NJ+2, NJP3=NJ+3, NJP4=NJ+4,  
     Q     NJM1=NJ-1, NJM2=NJ-2, NJM3=NJ-3, NJM4=NJ-4,  
     Q     NJQ2=2*NJ, NJQ2P1=NJQ2+1, NJQ2P2=NJQ2+2, NJ2=NJ/2)     
      PARAmeter (
     Q     NKP1=NK+1, NKP2=NK+2, NKP3=NK+3, NKP4=NK+4,  
     Q     NKM1=NK-1, NKM2=NK-2, NKM3=NK-3, NKM4=NK-4,  
     Q     NK2=NK/2, NK2P1=NK2+1, nk2p2=nk2+2,
     Q     nk2p3=nk2+3,  nk4=nk/4, nk8=nk/8) 
      PARAmeter (  
     Q     liP1=li+1, liP2=li+2, liP3=li+3, liP4=li+4,    
     Q     liM1=li-1, liM2=li-2, liM3=li-3, liM4=li-4,  
     Q     liQQ=li+20, lip5=li+5)                       
      PARAmeter (  
     Q     ljP1=lj+1, ljP2=lj+2, ljP3=lj+3, ljP4=lj+4,    
     Q     ljM1=lj-1, ljM2=lj-2, ljM3=lj-3, ljM4=lj-4,  
     Q     ljQQ=lj+20, ljp5=lj+5)                       
      PARAmeter (NOUT=ni,  NOUTP1=NOUT+1)                   
      parameter (nx=njp3, ny=nk2p3, kk=9, nxm2=nx-2, nym2=ny-2)
      PARAmeter (NTOT=NI*NJ*NK,                            
     Q     NSCALAR = NIP1*NJP1*NKP1,
     Q     NEDGEI=NI*NJP1*NKP1,NEDGEJ=NIP1*NJ*NKP1,NEDGEK=NIP1*NJP1*NK,
     Q     NFACEI=NIP1*NJ*NK,NFACEJ=NI*NJP1*NK,NFACEK=NI*NJ*NKP1,
     Q     NVAR0=11*NTOT + NI*NJ + NJ*NK + NI*NK,
     Q     NVAR1=11*NTOT + NI*NJ + NJ*NK + NI*NK+4*nX*NY ,
     $	   nvartot = nvar1 +2*njp1*nk2p1 +nedgei+nedgej+nedgek,
     Q     NVARALL = 2+ NVARTOT,
     Q     nVARBASE=8*ntot,
     Q     BIG2D=NI*NJ,
     Q     NVAROUT = 8*NTOT, NCURPOTOUT=4*NJP3*NK2p3, 
     Q     NPSIOUT = 2*NJP1*NKP1, NCONDOUT=4*NJP3*NK2P3,
     Q	   NIONVAR = 12*NJP1*NK2P3)
      real rearth,rion,rionm
      parameter( rearth=6.38e8,rion=6.5e8, rionm=6.5e6 )
      real erg2ev,effmass,boltzman
      parameter( erg2ev=1.6022e-12,effmass=2.24e-24,boltzman=4.8e-9)
      



