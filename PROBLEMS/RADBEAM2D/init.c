
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


rho=RHOAMB;
T=TAMB;
E=calc_LTE_EfromT(T);
uint=calc_PEQ_ufromTrho(T,rho);

#ifdef BLOB
rho=RHOAMB*(1.+BLOBP*exp(-((xx-BLOBX)*(xx-BLOBX)+(yy)*(yy)+(zz-BLOBZ)*(zz-BLOBZ))/BLOBW/BLOBW));
//T=TAMB*RHOAMB/rho;
E=calc_LTE_EfromT(T);
uint=calc_PEQ_ufromTrho(T,rho);

#endif

ldouble V=0.;

//zaczynam jednak od profilu analitycznego:   
ldouble r=xx;
ldouble mD=PAR_D/(r*r*sqrt(2./r*(1.-2./r)));
ldouble mE=PAR_E/(powl(r*r*sqrt(2./r),GAMMA)*powl(1.-2./r,(GAMMA+1.)/4.));
V=sqrt(2./r)*(1.-2./r)           ;
ldouble W=1./sqrt(1.-V*V*gg[1][1]);
rho=PAR_D/(r*r*sqrt(2./r));
uint=mE/W;

#ifdef FLATBACKGROUND
rho=RHOAMB;
T=TAMB;
uint=calc_PEQ_ufromTrho(T,rho);
V=0.;
#endif

E=calc_LTE_Efromurho(uint,rho);

//test!
//V=0.7;

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

prad_zamo2ff(pp,pp,gg,GG,eup);
prad_ff2lab(pp,pp,&geom);

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
