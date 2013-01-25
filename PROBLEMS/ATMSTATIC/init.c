
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

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);
ldouble gg[4][5],GG[4][5],tup[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
calc_tetrades(gg,tup,tlo);

ldouble pp[NV],T;

/************************/


//printf("%Le %Le %e %e\n",xx,f,KAPPAES,FLUXLEFT); getchar();


ldouble rho0=1.;
ldouble r0=2.e6;
ldouble u0=0.0001;
ldouble r=xx;

rho=rho0*r0/r;
uint=u0*r0*r0/r/r;

E=exp(1.);

uint=(4*r*u0 - 4*r*u0*GAMMA - 2*r*rho0 + 2*r0*rho0 + r*r0*rho0*Log(-2 + r) - r*r0*rho0*Log(r) - r*r0*rho0*Log(-2 + r0) + r*r0*rho0*Log(r0))/(4*r - 4*r*GAMMA);

//printf("%e %e %e %e\n",r,rho,E,uint);  getchar();

pp[0]=rho;
pp[1]=uint;
pp[2]=0.;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(rho,uint);
    
p2u(pp,uu,gg,GG);	 


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
