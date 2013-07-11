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
v6=fabs(get_cflag(ENTROPYFLAG,ix,iy,iz))+1.e-10;
v7=GAMMAM1*pp[UU];

