//user output - called from fileop.c

//v1 - 24
//..
//v7 - 30

v1=v2=v3=v4=v5=v6=v7=0.;

//tau in radius
ldouble tautot[3];
calc_tautot(pp,xxvec,dx,tautot);
v2=tautot[0];

//prad/ptotal
v5=get_cflag(HDFIXUPFLAG,ix,iy,iz);
v6=get_cflag(RADFIXUPFLAG,ix,iy,iz);
