//v1 - 24
//..
//v7 - 30

ldouble bcond[4],bcovd[4],bsqd;


#ifdef MAGNFIELD
calc_bcon_prim(pp,bcond,&geom);
indices_21(bcond,bcovd,gg); 
bsqd = dot(bcond,bcovd);
v1=bsqd/2.;
v2=calc_divB(ix,iy,iz);
#endif

#ifdef RADIATION
ldouble Rtt,Ehat,ucon[4],prad;
calc_ff_Rtt(pp,&Rtt,ucon,&geom);
Ehat=-Rtt; 
prad=Ehat/3.;
v3=prad;
#endif

v4=get_cflag(RADSOURCETYPEFLAG,ix,iy,iz);
