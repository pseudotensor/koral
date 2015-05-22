//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,
//	ldouble *uu,ldouble *pp,int ifinit,int BCtype)

/**********************/
//geometries
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

/**********************/

//radius
if(BCtype==XBCHI) //outflow in magn, atm in rad., atm. in HD
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
    if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
    if(ucon[1]<0.) //inflow, resseting to atmosphere
      {
	//atmosphere in rho,uint and velocities and zero magn. field
	//set_hdatmosphere(pp,xxvec,gg,GG,4);
	ucon[1]=0.;
	#ifdef MAGNFIELD
	pp[B2]=pp[B3]=0.;
	#endif
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	pp[VX]=ucon[1];
	pp[VY]=ucon[2];
	pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }

#ifdef RADIATION
    ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
    conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
    if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,urfcon,urfcon,MYCOORDS,BLCOORDS);
    if(urfcon[1]<0.) //inflow, resseting to atmosphere
      {
	//atmosphere in radiation
	//set_radatmosphere(pp,xxvec,gg,GG,0);
	urfcon[1]=0.;
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,urfcon,urfcon,BLCOORDS,MYCOORDS);
	conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
	pp[FX0]=urfcon[1];
	pp[FY0]=urfcon[2];
	pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
#endif

   
    

    p2u(pp,uu,&geom);
    return 0;  
  }
 else if(BCtype==XBCLO) //outflow near BH / black cylinder
   {
     iix=0;
     iiy=iy;
     iiz=iz;

      for(iv=0;iv<NV;iv++)
	{
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}

      if(RMIN>rhorizonBL) //do not allow for inflow
	{
	  //checking for the gas inflow
	  ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};    
	  conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	  
	  //if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
	  if(ucon[1]>0.) //inflow, resseting to atmosphere
	    {
	      //atmosphere in rho,uint and velocities and zero magn. field
	      //set_hdatmosphere(pp,xxvec,gg,GG,4);
	      ucon[1]=0.;
#ifdef MAGNFIELD
	      pp[B2]=pp[B3]=0.;
#endif
	      //  if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
	      conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	      pp[VX]=ucon[1];
	      pp[VY]=ucon[2];
	      pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
	    }

#ifdef RADIATION
	  ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
	  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
	  //if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,urfcon,urfcon,MYCOORDS,BLCOORDS);
	  if(urfcon[1]>0.) //inflow, resseting to atmosphere
	    {
	      //atmosphere in radiation
	      //set_radatmosphere(pp,xxvec,gg,GG,0);
	      urfcon[1]=0.;
	      //if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,urfcon,urfcon,BLCOORDS,MYCOORDS);
	      conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
	      pp[FX0]=urfcon[1];
	      pp[FY0]=urfcon[2];
	      pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
	    }
#endif

	}
      
     p2u(pp,uu,&geom);
     return 0;
   }

//reflections/outflow in theta 
//in 3D polar cells overwritten with #define CORRECT_POLARAXIS_3D
if(BCtype==YBCLO) //upper spin axis 
  {      
    
    if(MYCOORDS!=CYLCOORDS) //reflection from the spin axis
      {
	iiy=-iy-1;
	iiz=iz;
	iix=ix;
	gdet_src=get_g(g,3,4,iix,iiy,iiz);  
	gdet_bc=get_g(g,3,4,ix,iy,iz);  
	for(iv=0;iv<NV;iv++)
	  {
	    //v_theta
	    if(iv==VY || iv==B2 || iv==FY0)
	      pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	    else
	      pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }
      }
    else //cylindrical, outflow through lower boundary
      {
	iix=ix;
	iiy=0;
	iiz=iz;
    
	//copying everything
	for(iv=0;iv<=NV;iv++)
	  {
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }

	//checking for the gas inflow
	ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};    
	conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
	if(ucon[2]>0.) //inflow, resseting to atmosphere
	  {
	    //atmosphere in rho,uint and velocities and zero magn. field
	    //set_hdatmosphere(pp,xxvec,gg,GG,4);
	    ucon[2]=0.;
	    	#ifdef MAGNFIELD
	pp[B1]=pp[B3]=0.;
	#endif
	    if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
	    conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	    pp[VX]=ucon[1];
	    pp[VY]=ucon[2];
	    pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
	  }

#ifdef RADIATION
	ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
	conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,urfcon,urfcon,MYCOORDS,BLCOORDS);
	if(urfcon[2]>0.) //inflow, resseting to atmosphere
	  {
	    //atmosphere in radiation
	    //set_radatmosphere(pp,xxvec,gg,GG,0);
	    urfcon[2]=0.;
	    if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,urfcon,urfcon,BLCOORDS,MYCOORDS);
	    conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
	    pp[FX0]=urfcon[1];
	    pp[FY0]=urfcon[2];
	    pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
	  }
#endif
      } //end of cylindrical



    p2u(pp,uu,&geom);
    return 0;
  }

if(BCtype==YBCHI) //lower spin axis
  {
    if(MYCOORDS!=CYLCOORDS)
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
 
      }
    else //cylindrical, outflow through lower boundary
      {
	iix=ix;
	iiy=NY-1;
	iiz=iz;
    
	//copying everything
	for(iv=0;iv<=NV;iv++)
	  {
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }

	//checking for the gas inflow
	ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};    
	conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
	if(ucon[2]<0.) //inflow, resseting to atmosphere
	  {
	    //atmosphere in rho,uint and velocities and zero magn. field
	    //set_hdatmosphere(pp,xxvec,gg,GG,4);
	    ucon[2]=0.;
	    if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
	    conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	    pp[VX]=ucon[1];
	    pp[VY]=ucon[2];
	    pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
	  }

#ifdef RADIATION
	ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
	conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,urfcon,urfcon,MYCOORDS,BLCOORDS);
	if(urfcon[2]<0.) //inflow, resseting to atmosphere
	  {
	    //atmosphere in radiation
	    //set_radatmosphere(pp,xxvec,gg,GG,0);
	    urfcon[2]=0.;
	     	#ifdef MAGNFIELD
	pp[B1]=pp[B3]=0.;
	#endif
	    if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,urfcon,urfcon,BLCOORDS,MYCOORDS);
	    conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
	    pp[FX0]=urfcon[1];
	    pp[FY0]=urfcon[2];
	    pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
	  }
#endif
      } //end of cylindrical


    p2u(pp,uu,&geom); 
    return 0; 
  }
   
//periodic in phi:
iiz=iz;
iiy=iy;
iix=ix;
if(BCtype==ZBCLO) iiz=iz+NZ;
if(BCtype==ZBCHI) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

