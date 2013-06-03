
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
ldouble gg[4][5],GG[4][5],eup[4][4],elo[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(tmulo,ix,iy,iz,tlo);


ldouble pp[NV],T;

/************************/

ldouble rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut;


//flat gas profiles
Tgas=T_AMB;
rho=RHO_AMB;
uint=calc_PEQ_ufromTrho(Tgas,rho);


Fz=Fy=Fx=0.;
pp[0]=rho;
pp[1]=uint;
pp[2]=0.;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(rho,uint);

#ifdef RADIATION

ldouble E1,E2;

E=EEAMB;

pp[6]=E;
pp[7]=Fx;
pp[8]=Fy;
pp[9]=Fz; 

int irf;
#ifdef MULTIRADFLUID

for(irf=1;irf<NRF;irf++)
  {
    pp[FX(irf)]=pp[FY(irf)]=pp[FZ(irf)]=0.;
    pp[EE(irf)]=SMALL;
  }

#endif



#endif

prad_ff2lab(pp,pp,&geom);
p2u(pp,uu,&geom);	 

//print_Nvector(pp,NV);
#ifdef MULTIRADFLUID

redistribute_radfluids(pp,uu,&geom);
//print_Nvector(uu,NV);

u2p_rad(uu,pp,&geom,&irf);
#endif
//print_Nvector(pp,NV);getchar();


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
