
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

ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecSPH[4];

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecSPH,MYCOORDS,SPHCOORDS);
xx=xxvecSPH[1];
yy=xxvecSPH[2];
zz=xxvecSPH[3];


ldouble gg[4][5],GG[4][5];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomSPH;
fill_geometry_arb(ix,iy,iz,&geomSPH,SPHCOORDS);

ldouble pp[NV],ppback[NV],T;

//working in SPH
ldouble ggSPH[4][5],GGSPH[4][5];
calc_g_arb(xxvecSPH,ggSPH,SPHCOORDS);
calc_G_arb(xxvecSPH,GGSPH,SPHCOORDS);
ldouble eupSPH[4][4],eloSPH[4][4];
ldouble tupSPH[4][4],tloSPH[4][4];
calc_tetrades(ggSPH,tupSPH,tloSPH,SPHCOORDS);
calc_ZAMOes(ggSPH,eupSPH,eloSPH,SPHCOORDS);


ldouble podpierd=-(GGSPH[0][0]-2.*ELL*GGSPH[0][3]+ELL*ELL*GGSPH[3][3]);
ldouble ut=-1./sqrt(podpierd);

ut/=UTPOT; //rescales rin
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;
if(1)
  {
    //ambient
    set_hdatmosphere(pp,xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(pp,xxvec,geom.gg,geom.GG,0);

#endif
  }

pp[5]=calc_Sfromu(pp[0],pp[1]);


p2u(pp,uu,&geom);


//calculate the intensities
double RijM1[4][4];double M1[5];
calc_Rij_M1(pp,&geom,RijM1);
//converting to RADCLOSURECOORDS
trans22_coco(geom.xxvec, RijM1, RijM1, MYCOORDS, RADCLOSURECOORDS);
//to ortonormal
trans22_cc2on(RijM1,RijM1,tupSPH);
//input
M1[0]=RijM1[0][0];
M1[1]=RijM1[0][1];
M1[2]=RijM1[0][2];
M1[3]=RijM1[0][3];
M1[4]=pp[EE0];
      
#ifdef myVET
ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
#endif
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


//print_primitives(pp);
