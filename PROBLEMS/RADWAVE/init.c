
/*
  int
  set_initial_profile()
  {
  int ix,iy,iz;
  for(iz=0;iz<NZ;iz++)
  {
  for(iy=0;iy<NY;iy++)
  {
  for(ix=0;ix<NX;ix++)
  {
*/

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE,vx;  
ldouble xx,yy,zz;
ldouble uu[NV];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);
ldouble gg[4][5],GG[4][5],tup[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

ldouble pp[NV],T;

/************************/

ldouble t=0.;

//Jiang+12 waves
#if (NWAVE==5)

//printf("RHOZERO = %Le\nUINT = %Le\nT = %Le\nERAD = %Le\nSIGMA = %Le\n",(ldouble)RHOZERO,(ldouble)UINT,(ldouble)TEMP,(ldouble)ERAD,(ldouble)SIGMA_RAD);getchar();


rho=RHOZERO*(1+DRRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DRIM/DRRE*sin(OMRE*t-KK*xx)));
//ldouble DURE=DPRE/(GAMMA-1.); ldouble DUIM=DPIM/(GAMMA-1.);
uint=UINT*(1.+DURE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DUIM/DURE*sin(OMRE*t-KK*xx))) ;
ldouble cs=1/CC;
vx=0.+DVRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DVIM/DVRE*sin(OMRE*t-KK*xx)) ; //DVRE absolute!
E=ERAD*(1+DERE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DEIM/DERE*sin(OMRE*t-KK*xx)));
Fx=0.+ERAD*DFRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DFIM/DFRE*sin(OMRE*t-KK*xx));
Fz=Fy=0.;

//rho=RHOZERO;
//uint=UINT;
//E=ERAD;
//vx=0;
//Fx=0.;
#endif

//hydro density wave
#if (NWAVE==1)
rho=RHOZERO*(1.+AAA*cos(KK*xx));
uint=UINT;
vx=VVX;
#endif

//radiative hydro density wave
#if (NWAVE==3)
rho=RHOZERO*(1.+AAA*cos(KK*xx));
uint=UINT;
vx=VVX;
E=ERAD;
Fx=Fz=Fy=0.;
#endif

//hydro sound wave
#if (NWAVE==2)
rho=RHOZERO*(1.+AAA*cos(KK*xx));
uint=UINT*(1.+GAMMA*AAA*cos(KK*xx));
ldouble cs=1./CC;
vx=AAA*cos(KK*xx)*cs;
E=ERAD;
Fx=Fz=Fy=0.;
#endif

//radiative sound wave
#if (NWAVE==4)
rho=RHOZERO*(1.+GASFACTOR*AAA*cos(KK*xx));
uint=UINT*(1.+GASFACTOR*GAMMA*AAA*cos(KK*xx));
ldouble cs=1./CC;
vx=GASFACTOR*AAA*cos(KK*xx)*cs;
E=ERAD*(1.+ERADFACTOR*AAA*cos(KK*xx));
Fx=Fz=Fy=0.;
#endif

pp[0]=rho;
pp[1]=uint;
pp[2]=vx;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
pp[6]=E;
pp[7]=Fx;
pp[8]=Fy;
pp[9]=Fz; 

prad_ff2lab(pp,pp,&geom);
#endif

//print_Nvector(pp,NV); getchar();
p2u(pp,uu,&geom);	 


/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);

//if(isnan(get_u(p,5,ix,iy,iz))) {printf("pr: %d %d %d S: %Le\n",ix,iy,iz,0.);getchar();}

//mark initialy succesfull u2p_hot step
set_cflag(0,ix,iy,iz,0);
