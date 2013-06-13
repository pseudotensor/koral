//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

int iix,iiy,iiz,iv;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
/***********************************************/

/***********************************************/
/***********************************************/
if(ix>=NX) //analytical solution within the torus and atmosphere outside
  {
    ldouble theta=geomBL.yy;
    if(theta<INJTHETA) // atmosphere
      {
	set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,4);
	//rad atmosphere
#ifdef RADIATION
	set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
 
	//testing if interpolated primitives make sense
	check_floors_hd(pp,VELPRIM,&geom);
	//end of floor section

	p2u(pp,uu,&geom);
      }
    else //injection region
      {
	ldouble R=geomBL.xx;						
	ldouble Rt=R/2.;
	ldouble vr=-7.6e8 * ALPHAHDVISC * (MDOT*16.)*(MDOT*16.) * sqrt(1./Rt/Rt/Rt/Rt/Rt) * (1.-sqrt(3/Rt));
	ldouble MdotEdd = 2.23e18*MASS; //g/s
	ldouble Mdot = MDOT * MdotEdd;
	ldouble Rcgs = lenGU2CGS(R);
	ldouble Sigma = -Mdot / (2.*M_PI*vr*Rcgs);	
	ldouble rho0 = 35./16. * Sigma / ((M_PI/2. - INJTHETA) * lenGU2CGS(ROUT)) / 2.; 
	ldouble thetat = 1. - (theta - INJTHETA) / (M_PI/2. - INJTHETA);
	if(thetat>1.) thetat=1.;
	ldouble rho = rho0 * pow(1. - thetat*thetat,3.);
	rho = rhoCGS2GU(rho);
	ldouble uint = 0.0001*rho;
	pp[0]=rho;
	pp[1]=uint;

	ldouble uphi=(RKEP*RKEP/(sqrt(RKEP*(RKEP*RKEP-3.*RKEP)))); 
	ldouble ucon[4]={0.,velCGS2GU(vr),0.,uphi/R/R};
	conv_vels(ucon,ucon,VEL3,VEL4,geomBL.gg,geomBL.GG);
	trans2_coco(geomBL.xxvec,ucon,ucon,KERRCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);

	
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];

#ifdef RADIATION
	pp[6]=1.*uint;
	pp[7]=0.;
	pp[8]=0.;
	pp[9]=0.;
	
	//transforming rad primitives from BL to MYCOORDS
	trans_prad_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,geomBL.gg,geomBL.GG,geom.gg,geom.GG);
#endif 

	p2u(pp,uu,&geom);
      }
    return 0;
  }
 
/***********************************************/
/***********************************************/
if(ix<0) //outflow near BH
   {
     //extrapolating along MYCOORDS
     //gc radialcoordinate
     ldouble r=get_x(ix,0);
     //iix=0
     ldouble r0=get_x(0,0);
     ldouble r1=get_x(1,0);     

     for(iv=0;iv<NV;iv++)
       {
	 //linear extrapolation
	 pp[iv]=get_u(p,iv,0,iy,iz)+(get_u(p,iv,1,iy,iz)-get_u(p,iv,0,iy,iz))*(r-r0)/(r1-r0);
	 //constant
	 //pp[iv]=get_u(p,iv,0,iiy,iiz);
       }
 
      //testing if interpolated primitives make sense
      check_floors_hd(pp,VELPRIM,&geom);
      //end of floor section

      p2u(pp,uu,&geom);
      return 0;
   }

/***********************************************/
/***********************************************/
//reflections/outflow in theta 
if(iy<0.) //spin axis 
  {      
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    for(iv=0;iv<NV;iv++)
      {
	//theta vector components
	if(iv==3 || iv==8)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
    
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section

    p2u(pp,uu,&geom);
    return 0;
  }

/***********************************************/
/***********************************************/
if(iy>=NY) //equatorial plane
  {
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;
  	  
    for(iv=0;iv<NV;iv++)
      {
	if(iv==3 || iv==8)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section

    p2u(pp,uu,&geom); 
    return 0; 
  }
   
/***********************************************/
/***********************************************/
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

//and that is all
 
return 0;

