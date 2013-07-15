//user output - called from fileop.c

//v1 - 24
//..
//v7 - 30

v1=v2=v3=v4=v5=v6=v7=0.;

//tau in radius
ldouble tautot[3];
calc_tautot(pp,xxvec,dx,tautot);
v2=tautot[1];

//prad/ptotal
v5=get_cflag(HDFIXUPFLAG,ix,iy,iz);
v6=get_cflag(RADFIXUPFLAG,ix,iy,iz);

ldouble shear[4][4],shearon[4][4];
ldouble vdiff2,nu;
calc_nu_shearviscosity(pp,&geom,shear,&nu,&vdiff2);
ldouble ev[4],evmax;
trans22_cc2on(shear,shearon,geom.tup);
evmax=calc_eigen_4x4(shearon,ev);
//xxvec[4]={0.,geom.xx,geom.yy,geom.zz};
ldouble xxvecBL[4];
coco_N(geom.xxvec,xxvecBL,MYCOORDS,BLCOORDS);
v1=evmax/pow(xxvecBL[1],-.5);

rho=pp[RHO];
ldouble uint2=pp[UU];
 
ldouble pre=(GAMMA-1.)*uint2;
ldouble cs2=GAMMA*pre/(rho+uint2+pre);

v2=sqrt(cs2);
