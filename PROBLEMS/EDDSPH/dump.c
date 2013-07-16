//user output - called from fileop.c

//v1 - 24
//..
//v7 - 30

ldouble Rij[4][4];
calc_Rij(pp,&geom,Rij);

trans22_cc2on(Rij,Rij,geom.tup);

v1=Rij[0][0];
v2=Rij[0][1];///v1;
v3=Rij[1][1];///v1;
v4=Rij[2][2];///v1;
v5=Tgas;
v6=Trad;

//tau in radius
ldouble tautot[3];
calc_tautot(pp,geom.xxvec,dx,tautot);

ldouble kappa=calc_kappa(pp[RHO],Tgas,-1.,-1.,-1.);
ldouble kappaes=calc_kappaes(pp[RHO],Tgas,-1.,-1.,-1.);

v3=kappa*dx[1];
v4=kappaes*dx[1];

