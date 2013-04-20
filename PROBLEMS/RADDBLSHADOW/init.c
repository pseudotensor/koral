
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

struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV];

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);
ldouble gg[4][5],GG[4][5],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(tmulo,ix,iy,iz,tlo);
ldouble pp[NV],T;

/*****************************/

ldouble pamb=calc_PEQ_ufromTrho(TAMB,RHOAMB);
#ifdef BLOB
rho=(RHOBLOB-RHOAMB)*exp(-powl(sqrtl(xx*xx+yy*yy+zz*zz),2.)/BLOBW)+RHOAMB;      
#else
rho=RHOAMB;
#endif
mx=my=mz=0.;
T=TAMB*RHOAMB/rho;
//T=TAMB;

uint=calc_PEQ_ufromTrho(T,rho);
E=calc_LTE_EfromT(T);

Fx=0.;
Fy=0.;
Fz=0.;	      	      

ldouble V=0.;


Fx=Fy=Fz=0.*E;
pp[0]=rho;
pp[1]=uint;
pp[2]=-V;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(rho,uint);
pp[6]=E;
pp[7]=Fx;
pp[8]=Fy;
pp[9]=Fz;

prad_ff2lab(pp,pp,&geom);

/*
pp[6]=1.;
pp[7]=pp[8]=0.;
*/

p2u(pp,uu,&geom);

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
