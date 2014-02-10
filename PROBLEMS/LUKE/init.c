
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
ldouble gg[4][5],eup[4][4],elo[4][4],GG[4][5];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

ldouble pp[NV],T;

/************************/

pp[0]=1.;
pp[1]=UINTFRAC;
pp[2]=pp[3]=pp[4]=0.;

#ifdef SPHERICAL
pp[2]=VELINX*cos(zz)*sin(yy)/sqrt(gg[1][1]);
pp[3]=VELINX*cos(zz)*cos(yy)/sqrt(gg[2][2]);
pp[4]=-VELINX*sin(zz)/sqrt(gg[3][3]);
#endif

#ifdef CYLINDRICAL
pp[2]=VELINX*cos(zz);
pp[4]=-VELINX*sin(zz)/sqrt(gg[3][3]);
#endif

pp[5]=calc_Sfromu(pp[0],pp[1]);	      

//converting from 3vel to VELPRIM
conv_velsinprims(pp,VEL3,VELPRIM,gg,GG);
	      
if(ix>=-1) //conserved required for ix=-1 only
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
