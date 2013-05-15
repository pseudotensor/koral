
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
ldouble uu[NV],xxvec[4],xxvecBL[4];

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,MYCOORDS);
xx=xxvecBL[1];
yy=xxvecBL[2];
zz=xxvecBL[3];


ldouble gg[4][5],GG[4][5],eup[4][4],elo[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,MYCOORDS);

pick_T(tmulo,ix,iy,iz,tlo);
calc_ZAMOes(gg,eup,elo,MYCOORDS);

ldouble pp[NV],ppback[NV],T;

//working in BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,MYCOORDS);
calc_G_arb(xxvecBL,GGBL,MYCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,MYCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,MYCOORDS);



ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;
if(1)
  {
    //ambient
    set_hdatmosphere(pp,xxvec,gg,GG,0);
    pp[2]=0.;
    pp[0]=1.;
    pp[1]=.1;
#ifdef RADIATION
    set_radatmosphere(pp,xxvec,gg,GG,0);
    pp[6]=0.0001;
#ifdef LABRADFLUXES
    pp[6]=-0.0001;
#endif
    pp[7]=pp[8]=pp[9]=0.;

    int irf;
    for(irf=1;irf<NRF;irf++)
      {
	pp[EE(irf)]=EEFLOOR;
	pp[FX(irf)]=0.;
	pp[FY(irf)]=0.;
	pp[FZ(irf)]=0.;
      }

    /*
    pp[6]=ERADATMMIN;
    pp[7]=0.;
    pp[8]=0.;
    pp[9]=0.;

    //transforming BL ZAMO radiative primitives to BL non-ortonormal primitives
    ldouble eupBL[4][4],eloBL[4][4];
    ldouble tupBL[4][4],tloBL[4][4];
    calc_tetrades(ggBL,tupBL,tloBL,MYCOORDS);
    calc_ZAMOes(ggBL,eupBL,eloBL,MYCOORDS);
    prad_zamo2ff(pp,pp,ggBL,GGBL,eupBL);
    prad_ff2lab(pp,pp,ggBL,GGBL,tloBL);
    //transforming radiative primitives from BL to MYCOORDS
    trans_prad_coco(pp, pp, MYCOORDS, MYCOORDS,xxvec,ggBL,GGBL,gg,GG);
    */
    
#endif
  }

pp[5]=calc_Sfromu(pp[0],pp[1]);

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section

p2u(pp,uu,&geom);


#ifdef MULTIRADFLUID

int irf;
redistribute_radfluids(pp,uu,&geom);
u2p_rad(uu,pp,&geom,&irf);

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
