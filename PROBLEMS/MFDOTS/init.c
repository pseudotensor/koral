
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

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);
ldouble gg[4][5],GG[4][5],tup[4][4],tlo[4][4],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(tmuup,ix,iy,iz,tup);
pick_T(tmulo,ix,iy,iz,tlo);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
ldouble pp[NV],T;

/*****************************/

pp[0]=1.;
pp[1]=1.;
pp[2]=0.;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(rho,uint);
pp[6]=LTEFACTOR*calc_LTE_Efromurho(pp[0],pp[1]);
pp[7]=0.;
pp[8]=0.;
pp[9]=0.;

int irf=0;

    if(ix==IXDOT1 && iy==IYDOT1 && iz==IZDOT1 && 1)
      {
	pp[6]*=1000.;
	pp[FX(irf)]=FXDOT1*pp[6];
	pp[FY(irf)]=FYDOT1*pp[6];
	pp[FZ(irf)]=FZDOT1*pp[6];
      }
   if(ix==IXDOT2 && iy==IYDOT2 && iz==IZDOT2 && 1)
      {
	pp[6]*=1000.;
	pp[FX(irf)]=FXDOT2*pp[6];
	pp[FY(irf)]=FYDOT2*pp[6];
	pp[FZ(irf)]=FZDOT2*pp[6];
      }

for(irf=1;irf<NRF;irf++)
  {
    pp[EE(irf)]=EEFLOOR;
    pp[FX(irf)]=0.;
    pp[FY(irf)]=0.;
    pp[FZ(irf)]=0.;
  }

prad_ff2lab(pp,pp,&geom);

p2u(pp,uu,&geom );
//printf("%d %d %d\n",ix,iy,iz);
//print_Nvector(uu,NV);
//print_Nvector(uu,NV);getchar();

#ifdef MULTIRADFLUID

//  print_Nvector(uu,NV);
redistribute_radfluids(pp,uu,&geom);
//    print_Nvector(uu,NV);
u2p_rad(uu,pp,&geom,&irf);

#endif

/**************************/

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
