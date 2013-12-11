int diskatboundary(ldouble *pp, void *ggg, void *gggBL);

/**********************/

//definitions
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecBL[4],xx,yy,zz;

//coordinates
get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

//metric
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);

//BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,KERRCOORDS);
calc_G_arb(xxvecBL,GGBL,KERRCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);
/**********************/

//outer edge, outflows with velocity check
if(ix>=NX) 
  {
    if(diskatboundary(pp, &geom, &geomBL)<0)
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
	if(ucon[1]<0.) //inflow, resetting the radial velocity
	  {
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
	if(urfcon[1]<0.) //inflow, resetting the radial velocity
	  {
	    //set_radatmosphere(pp,xxvec,gg,GG,0);
	    urfcon[1]=0.;
	    conv_vels(urfcon,urfcon,VEL4,VELPRIM,geom.gg,geom.GG);
	    pp[FX0]=urfcon[1];
	    pp[FY0]=urfcon[2];
	    pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
	  }
#endif
      }
    p2u(pp,uu,&geom);
	
    return 0;  
  }

 else if(ix<0) //outflow inside BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     //linear extrapolation
      for(iv=0;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       }

     p2u(pp,uu,&geom);
     return 0;
   }

//reflections in theta 
if(iy<0.) //axis
  {      
    
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	//theta component
	if(iv==VY || iv==B2 || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }

    p2u(pp,uu,&geom);
    return 0;
  }
if(iy>=NY) //axis or eq.plane
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

return 0;

