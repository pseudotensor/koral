//KORAL - init.c
//sets initial state
//**********************
//included within:
/*
  int set_initial_profile()
  {
  int ix,iy,iz;
  for(iz=0;iz<NZ;iz++)
  {
  for(iy=0;iy<NY;iy++)
  {
  for(ix=0;ix<NX;ix++)
  {
*/
//from problem.c
//**********************


  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV];
ldouble gg[4][5],GG[4][5],tup[4][4],tlo[4][4];
ldouble pp[NV],T;

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);

pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
calc_tetrades(gg,tup,tlo,MYCOORDS);

/************************/
/************************/
/************************/
/* modify below */
      
ldouble rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut,ux;




if(xx<(MAXX+MINX)/2.)
  {
    rho=1.;
    if(NTUBE==1) {uint = 3.e-5 / (GAMMA - 1.); E=1.e-8; Fx=1.e-2*E;ux=0.015;}
    if(NTUBE==2) {uint = 4.e-3 / (GAMMA - 1.);E=2.e-5; Fx=1.e-2*E;ux=0.25;}
    if(NTUBE==3 || NTUBE==31) {uint = 60. / (GAMMA - 1.);E=2.; Fx=1.e-2*E;ux=10.;}
    if(NTUBE==4 || NTUBE==41) {uint = 6.e-3 / (GAMMA - 1.);E=0.18; Fx=1.e-2*E;ux=0.69;}	  
    if(NTUBE==5) {uint = 60. / (GAMMA - 1.);E=2.; Fx=1.e-2*E;ux=1.25;}
  }
 else
   {
     if(NTUBE==1) {rho=2.4;uint = 1.61e-4/ (GAMMA - 1.); E=2.51e-7; Fx=1.e-2*E;ux=6.25e-3;}
     if(NTUBE==2) {rho=3.11;uint = 0.04512 / (GAMMA - 1.);E=3.46e-3; Fx=1.e-2*E;ux=0.0804;}
     if(NTUBE==3 || NTUBE==31) {rho=8.0;uint = 2.34e3 / (GAMMA - 1.);E=1.14e3; Fx=1.e-2*E;ux=1.25;}
     if(NTUBE==4 || NTUBE==41) {rho=3.65;uint =3.59e-2 / (GAMMA - 1.);E=1.30; Fx=1.e-2*E;ux=0.189;}	  
     if(NTUBE==5) {rho=1.0;uint = 60. / (GAMMA - 1.);E=2.; Fx=1.e-2*E;ux=1.10;}
   }



ut=sqrt(ux*ux+1.);
vx=ux/ut;
Fz=Fy=0.;
pp[0]=rho;
pp[1]=uint;

    //pp[2]=vx;
    //3-velocity in BL transformed to VELR
    ldouble ucon[4]={0.,vx,0.,0.};
    conv_vels(ucon,ucon,VEL3,VELPRIM,gg,GG);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3]; 

pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
pp[6]=E;
pp[7]=Fx;
pp[8]=Fy;
pp[9]=Fz; 

prad_ff2lab(pp,pp,&geom);

//pp[2]=0.;
#endif

//converting to conserved
p2u(pp,uu,&geom);	

/*
printf("%d %d %d\n",ix,iy,iz);
print_Nvector(uu,NV);

int t1[2],t2[3];
u2p(uu,pp,&geom,t1,t2);
print_Nvector(pp,NV);

p2u(pp,uu,&geom);	
getchar();
*/
/* modify above */
/***********************************************/
/***********************************************/
/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
//mark initialy succesfull u2p_hot step
set_cflag(0,ix,iy,iz,0);
