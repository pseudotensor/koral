
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
ldouble gg[4][5],GG[4][5],tup[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
ldouble pp[NV],T;

/************************/

ldouble rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut;

//at outern boundary
ldouble f = (ldouble)KAPPAES*FLUXLEFT*MINX*MINX;


//printf("flux: %e\n",fluxGU2CGS(FLUXLEFT)*4.*Pi*GMC2CM*GMC2CM*MINX*MINX);getchar();


ldouble p0=K_BOLTZ*RHOAMB*TAMB/MU_GAS/M_PROTON;	      
ldouble KKK=p0/powl(RHOAMB,GAMMA);
ldouble C3=GAMMA*KKK/(GAMMA-1.)*powl(RHOAMB,GAMMA-1.)-(1.-f)*(1./MINX+0.*1./MINX/MINX+0.*4./3./MINX/MINX/MINX);

pp[0]=powl(GAMMAM1/GAMMA/KKK*(C3+(1.-f)*(1./xx+0.*1./xx/xx+0.*4./3./xx/xx/xx)),1./GAMMAM1);

ldouble pre=KKK*powl(pp[0],GAMMA);

pp[1]=pre/GAMMAM1;

Fz=Fy=0.;
Fx=FLUXLEFT*(MINX/xx)*(MINX/xx);

#ifdef THIN
E=Fx/FERATIO;
#endif

#ifdef THICK
E=calc_LTE_EfromT(calc_PEQ_Tfromrhou(pp[0],pp[1]));
#endif




//printf("%Le %Le %e %e\n",xx,f,KAPPAES,FLUXLEFT); getchar();
pp[2]=0.;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
pp[6]=E;
pp[7]=Fx;
pp[8]=Fy;
pp[9]=Fz; 

prad_ff2lab(pp,pp,&geom);

/*
printf("%d\n",ix);
print_Nvector(pp,NV);
getchar();
*/

#endif	      
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
