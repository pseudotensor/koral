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

/**********************/

//outer edge, outflows with velocity check
if(BCtype==XBCHI)
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
	ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]},uconBL[4];    
	conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	trans2_coco(geom.xxvec,ucon,uconBL,MYCOORDS,BLCOORDS);
	
     	/*checks if boyer-lindquist 1 coordinate (theta) is within the appropriate bound
	  for inflow of gas; if so, sets up the inflow*/

	if(uconBL[1]<0.) //inflow
	  {
	    //set_hdatmosphere(pp,xxvec,gg,GG,4); //atmosphere with zero BL velocity	    
	    uconBL[1]=0.;
	    trans2_coco(geomBL.xxvec,uconBL,ucon,BLCOORDS,MYCOORDS);
	    conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	    pp[VX]=ucon[1];
	    pp[VY]=ucon[2];
	    pp[VZ]=ucon[3];
	  }

#ifdef RADIATION
	ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]},urfconBL[4];    
	conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
	trans2_coco(geom.xxvec,urfcon,urfconBL,MYCOORDS,BLCOORDS);
	if(urfconBL[1]<0.) //inflow, resetting the radial velocity
	  {
	    //set_radatmosphere(pp,xxvec,gg,GG,1);
	    urfconBL[1]=0.;
	    trans2_coco(geomBL.xxvec,urfconBL,urfcon,BLCOORDS,MYCOORDS);
	    conv_vels(urfcon,urfcon,VEL4,VELPRIM,geom.gg,geom.GG);
	    pp[FX0]=urfcon[1];
	    pp[FY0]=urfcon[2];
	    pp[FZ0]=urfcon[3];	    
	  }
#endif
      
#ifdef MAGNFIELD
    pp[B1]=pp[B2]=pp[B3]=0.;
#endif
      }
  
    p2u(pp,uu,&geom);

    /*
      print_primitives(pp);
      print_conserved(uu);getchar();
    */

    return 0;  
  }

 else if(BCtype==XBCLO) //outflow inside BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     //linear extrapolation
      for(iv=0;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       }

      if(RMIN>2.) //boundary outside BH
	{
	  //checking for the gas outflow
	  ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]},uconBL[4];    
	  conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	  trans2_coco(geom.xxvec,ucon,uconBL,MYCOORDS,BLCOORDS);
	
	  if(uconBL[1]>0.) //outflow
	    {
	      uconBL[1]=0.;
	      trans2_coco(geomBL.xxvec,uconBL,ucon,BLCOORDS,MYCOORDS);
	      conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	      pp[VX]=ucon[1];
	      pp[VY]=ucon[2];
	      pp[VZ]=ucon[3];
	    }
	}

     p2u(pp,uu,&geom);
     return 0;
   }

//reflections in theta 
if(BCtype==YBCLO) //axis
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
if(BCtype==YBCHI)//axis or eq.plane
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

