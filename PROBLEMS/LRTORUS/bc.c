//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

/**********************/
//geometries
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecBL[4],xx,yy,zz;

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);

//working in BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,KERRCOORDS);
calc_G_arb(xxvecBL,GGBL,KERRCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);
/**********************/

//radius
if(ix>=NX) //outflow in magn, atm in rad., atm. in HD
  {
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
    //copying everything
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    //checking for the gas inflow
    ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};    
    conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
    if(ucon[1]<0.) //inflow, resseting to atmosphere
      {
	//atmosphere in rho,uint and velocities and zero magn. field
	//set_hdatmosphere(pp,xxvec,gg,GG,4);
	ucon[1]=0.;
	conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	pp[VX]=ucon[1];
	pp[VY]=ucon[2];
	pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }

#ifdef RADIATION
    ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
    conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
    if(urfcon[1]<0.) //inflow, resseting to atmosphere
      {
	//atmosphere in radiation
	//set_radatmosphere(pp,xxvec,gg,GG,0);
	urfcon[1]=0.;
	conv_vels(urfcon,urfcon,VEL4,VELPRIM,geom.gg,geom.GG);
	pp[FX0]=urfcon[1];
	pp[FY0]=urfcon[2];
	pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
#endif

    /*
    #ifdef RADIATION
    pp[EE0]=pp[UU]/1.e6;
    pp[FX0]=pp[VX];
    pp[FY0]=pp[VY];
    pp[FZ0]=pp[VZ];
    #endif
    */
    

    p2u(pp,uu,&geom);
    return 0;  
  }
 else if(ix<0) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     //gc
     ldouble r=xxvec[1];

     //iix=0
     ldouble xxout[4]={0.,get_x(0,0),get_x(iiy,1),get_x(iiz,2)};
     ldouble r0=xxout[1];   
     
     //iix=1
     xxout[1]=get_x(1,0);
     ldouble r1=xxout[1];   
     

     //linear extrapolation
      for(iv=0;iv<NV;iv++)
       {
	 //pp[iv]=get_u(p,iv,0,iiy,iiz)+(get_u(p,iv,1,iiy,iiz)-get_u(p,iv,0,iiy,iiz))*(r-r0)/(r1-r0);
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       }

      //      pp[0] = get_u(p,0,0,iiy,iiz)*pow(r/r0,-1.5);
      //pp[1] = get_u(p,1,0,iiy,iiz)*pow(r/r0,-1.5);
 

     p2u(pp,uu,&geom);
     return 0;
   }

//reflections/outflow in theta 
if(iy<0.) //spin axis 
  {      
    
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	//v_theta
#ifndef PUREAXISOUTFLOW
	if(iv==VY || iv==B2 || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
#endif
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
     

    p2u(pp,uu,&geom);
    return 0;
  }
if(iy>=NY) //equatorial plane
  {
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
  	  
    for(iv=0;iv<NV;iv++)
      {
	if(iv==VY || iv==B2 || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	
      }
 
    p2u(pp,uu,&geom); 
    return 0; 
  }
   
//periodic in phi:
iiz=iz;
iiy=iy;
iix=ix;
if(iz<0) iiz=iz+NZ;
if(iz>NZ-1) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

